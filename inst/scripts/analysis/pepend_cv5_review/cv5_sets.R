## The three sets the PCA facets use, built under ONE cleanliness rule.
##
## Gene-level exclusion ("this gene supplied a negative anywhere") is too
## blunt: GDF5 got one negative at C-anchor 348, whose window spans 321-356,
## while its own end windows are 375-410 (N) and 474-509 (C). The model never
## saw a window containing either true end. So the test is WINDOW OVERLAP --
## an end is unseen iff no training-negative window overlaps its own 36-mer.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]
source("R/pepend_windows.R")
A2 <- path.expand("~/AF2_analysis")
LIM <- 60                                   # peptide := cleavage product < 60 aa
W   <- LF_PEPEND$seq_len                    # 36

ll   <- as.data.frame(readRDS("inst/extdata/ligand_list.rds"))
x    <- readRDS(file.path(A2, "lf_pepend_nn_input_near.rds"))
kn   <- as.data.frame(x$knowns)
prec <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))$prec
im   <- readRDS("data/id_mapping.rds")
a2g  <- setNames(im[["Gene Names (primary)"]], im$Entry)
ll$gene <- unname(a2g[ll$accession]); ll$len <- ll$end - ll$start + 1

## every training-negative window, as a precursor interval
negs <- bind_rows(lapply(c("N", "C"), function(tm)
  bind_rows(lapply(c("train", "val"), function(sp) {
    d <- x$nn_input[[tm]][[sp]]; i <- d$known == 0
    if (!any(i)) return(NULL)
    s <- lf_pepend_start(tm, d$anchor[i])
    tibble(accession = d$accession[i], lo = s, hi = s + W - 1L)
  }))))
negs <- as.data.table(negs); setkey(negs, accession)
cat(sprintf("training-negative windows: %d over %d precursors\n",
            nrow(negs), uniqueN(negs$accession)))

## does any negative window overlap this end's own window?
overlapped <- function(acc, term, anchor) {
  lo <- lf_pepend_start(term, anchor); hi <- lo + W - 1L
  n <- negs[.(acc), nomatch = 0L]
  if (!nrow(n)) return(FALSE)
  any(n$lo <= hi & lo <= n$hi)
}

p <- ll %>% filter(!is.na(start), !is.na(end)) %>%
  distinct(accession, start, end, .keep_all = TRUE) %>%
  filter(accession %in% prec$accession, !gene %in% kn$gene)
j <- match(p$accession, prec$accession)
p$n_prot <- prec$n_prot[j]; p$c_prot <- prec$c_prot[j]
p$name <- ifelse(is.na(p$ligand_name) | !nzchar(p$ligand_name), p$gene, p$ligand_name)

## one row per END, so each terminus is judged on its own window
ends <- bind_rows(
  p %>% transmute(gene, accession, name, len, term = "N", anchor = start, n_prot, c_prot),
  p %>% transmute(gene, accession, name, len, term = "C", anchor = end,   n_prot, c_prot)) %>%
  filter(anchor >= n_prot, anchor <= c_prot)            # inside the scanned range
ends$dirty <- mapply(overlapped, ends$accession, ends$term, ends$anchor)
ends$set <- ifelse(ends$len < LIM, "never seen (peptide, < 60 aa)",
                                   "large cleavage product ends")
clean <- ends %>% filter(!dirty)

cat("\n=== ends available, by set ===\n")
print(clean %>% count(set, term) %>% tidyr::pivot_wider(names_from = term, values_from = n) %>%
      as.data.frame(), row.names = FALSE)
print(clean %>% group_by(set) %>% summarise(ends = n(), peptides = n_distinct(paste(accession, len)),
      genes = n_distinct(gene), med_len = median(len), .groups = "drop") %>% as.data.frame(),
      row.names = FALSE)
cat(sprintf("\ndropped for window overlap: %d of %d ends (%d precursors rescued that a\n",
            sum(ends$dirty), nrow(ends),
            length(setdiff(unique(ends$accession[!ends$dirty]),
                           unique(ends$accession[!ends$accession %in% negs$accession])))))
cat("gene-level rule would have thrown out)\n")
cat("\nGDF5 / GDF6:\n")
print(ends %>% filter(gene %in% c("GDF5","GDF6")) %>%
      transmute(gene, term, anchor, len, dirty, set) %>% as.data.frame(), row.names = FALSE)
write.csv(clean, file.path(S2, "cv5_unseen_ends.csv"), row.names = FALSE)
cat("\nwrote cv5_unseen_ends.csv\n")
