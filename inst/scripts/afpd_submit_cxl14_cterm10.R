#!/usr/bin/env Rscript
## CXL14,100-109 against the GPCR panel, in three priority tiers.
##
## The ligand is NAWNEKRRVY, the C-terminal 10-mer of CXCL14.  It is a window
## of a real UniProt entry, so it needs no features of its own: AlphaPulldown
## slices the range out of the full-length CXL14 pickle at prediction time.
##
## Numbering is UniProt, not mature.  CXCL14's signal peptide is 1-34, so the
## mass-spec annotation "66-75" for this fragment is mature numbering and
## would address a completely different region if used here.  The containing
## proform KLQSTKRFIKWYNAWNEKRRVYEE is CXL14,88-111, which is the same start
## used by afpd_submit_CXCL14_deeper.R.
##
## Tiers are the requested order, which is not the order in
## `ecb: Order of runs (priority)` -- GPR15 and GPR25 are #1 there but asked
## for second here.  Tier 1 and 2 get one receptor per job file so they run in
## parallel rather than sequentially inside one allocation; tier 3 is batched
## for throughput.  Job files are numbered in priority order, so --first/--last
## on afpd_submit_jobs_cli.R is a priority window.

.libPaths("/home/groups/ebutcher/programs/pipeline/R_libs4.1")
suppressPackageStartupMessages({
  library(ligandFinder); library(dplyr); library(purrr); library(tibble)
})

LIGAND   <- "CXL14,100-109"
RUN      <- "cxl14_cterm10"
JOB_DIR  <- file.path("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts", RUN)
OUT_DIR  <- file.path("/scratch/groups/ebutcher/deorphan/models", RUN)
GROUP_SIZE_BULK <- 48

TIER1 <- c("AGTR2", "AGTR1", "BKRB1")                        # AGTR2, AGTR1, BDKRB1
TIER2 <- c("RL3R1", "RL3R2", "APJ", "GPR15", "GPR25")        # RXFP3, RXFP4, APLNR, GPR15, GPR25

## --- the standard receptor set -------------------------------------------
gpcr_sub <- gpcr_list %>%
  filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  filter(`ecb: Prioritization Notes` != "Small organic molecule") %>%
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))

message("standard receptor set: ", nrow(gpcr_sub))

named <- c(TIER1, TIER2)
absent <- setdiff(named, gpcr_sub[["uniprot_name"]])
if (length(absent) > 0) {
  stop("requested receptor(s) not in the standard set: ", paste(absent, collapse = ", "),
       "\nThey are either not priority-flagged, have no BW table, or are small-molecule ",
       "receptors. Add them explicitly if that is intended.")
}

## --- tier assignment ------------------------------------------------------
receptors <- gpcr_sub %>%
  transmute(uniprot_name,
            receptor = model,
            tier = case_when(uniprot_name %in% TIER1 ~ 1L,
                             uniprot_name %in% TIER2 ~ 2L,
                             TRUE ~ 3L)) %>%
  mutate(tier_rank = case_when(tier == 1L ~ match(uniprot_name, TIER1),
                               tier == 2L ~ match(uniprot_name, TIER2),
                               TRUE ~ NA_integer_)) %>%
  arrange(tier, tier_rank, uniprot_name)

to_run <- receptors %>%
  mutate(ligand = LIGAND,
         model  = paste(receptor, ligand, sep = ";"))

message("tier sizes: ", paste(sprintf("tier %d = %d", 1:3,
        c(sum(to_run$tier == 1), sum(to_run$tier == 2), sum(to_run$tier == 3))),
        collapse = ", "))

## --- every token must resolve to a pickle, checked before sbatch ---------
chk <- afpd_check_features(afpd_job_names(to_run[["model"]]))
if (any(!chk[["present"]])) {
  stop("no features for: ",
       paste(unique(chk[["name"]][!chk[["present"]]]), collapse = ", "))
}

## --- job files, numbered in priority order ------------------------------
## one receptor per file for tiers 1-2 so they occupy separate allocations
to_run <- to_run %>%
  mutate(block = if_else(tier <= 2L, row_number(), NA_integer_))

n_fast <- sum(to_run[["tier"]] <= 2L)
bulk_idx <- seq_len(nrow(to_run) - n_fast)
to_run[["block"]][to_run[["tier"]] == 3L] <-
  n_fast + ceiling(bulk_idx / GROUP_SIZE_BULK)

to_run <- to_run %>% mutate(group = paste0("job", block, ".txt"))

dir.create(JOB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = file.path(JOB_DIR, .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))

write.table(to_run %>% select(uniprot_name, receptor, ligand, tier, group, model),
            file = file.path(JOB_DIR, "manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

n_jobs <- length(unique(to_run[["group"]]))
message("\nwrote ", n_jobs, " job file(s) to ", JOB_DIR)
message("output -> ", OUT_DIR, " (on scratch, which is where the metric pipeline looks)")

t1 <- range(to_run$block[to_run$tier == 1L])
t2 <- range(to_run$block[to_run$tier == 2L])
t3 <- range(to_run$block[to_run$tier == 3L])

cat(sprintf("
tier 1 (AGTR2, AGTR1, BDKRB1)         job%d-%d
tier 2 (RXFP3, RXFP4, APLNR, GPR15/25) job%d-%d
tier 3 (remaining %d receptors)        job%d-%d

Submit tier 1 now:
  Rscript $R_LIBS/ligandFinder/exec/afpd_submit_jobs_cli.R \\
    -i %s -o %s -f %d -l %d -m 16

then tier 2 (-f %d -l %d), then tier 3 (-f %d -l %d).
", t1[1], t1[2], t2[1], t2[2], sum(to_run$tier == 3L), t3[1], t3[2],
   JOB_DIR, OUT_DIR, t1[1], t1[2], t2[1], t2[2], t3[1], t3[2]))
