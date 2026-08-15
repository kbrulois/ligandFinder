## =============================================================================
## Build an initial-guess input: an OPEN receptor conformation with the test
## peptide docked into its pocket.
##
## Instead of the apo AF-DB monomer -- whose N-terminus and/or ECL2 are usually
## sitting in the orthosteric pocket -- we take a model where some peptide
## already got into the pocket (selected by ap_initial_guess_active_state_receptors.R,
## extracted by ap_extract_open_models.R). That model's receptor chain is an open,
## holo-like conformation.
##
## The RESIDENT peptide in that model is used only as a POSE TEMPLATE. The
## coordinates that actually go into the output come from a freshly cleaved
## AF-DB structure of the test peptide, superposed onto the resident peptide's
## pocket-binding terminus. So the pose is borrowed; the geometry is the test
## peptide's own.
##
## This replaces the geometric placement that position_ligand_initial_guess()
## did -- a real in-pocket pose beats a computed one. That function is still the
## fallback for receptors with no open model at all (see the coverage table from
## ap_initial_guess_active_state_receptors.R).
## =============================================================================

.libPaths("/home/groups/ebutcher/programs/pipeline/R_libs4.1")
library(ligandFinder)
suppressPackageStartupMessages(library(dplyr))

## ---- options ---------------------------------------------------------------
af_db        <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/af_db"
open_dir     <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/open_models"
manifest_csv <- file.path(open_dir, "open_receptors_manifest.csv")
out_dir      <- "/scratch/groups/ebutcher/kevin/ig_test"

## receptor to build for (matches `receptor` / `receptor_id` in the manifest)
rec_name   <- "NK2R"
rec_id     <- "P21452"

## test peptide: Neuropeptide K = TKN1 residues 72-107
lig_id     <- "P20366"
lig_start  <- 72
lig_end    <- 107
lig_seq    <- "DADSSIEKQVALLKALYGHGQISHKRHKTDSFVGLM"   # sanity check; NULL to skip
lig_label  <- "TKN1"

terminus   <- "C"    # which peptide terminus goes into the pocket (NKA family: C)
n_anchor   <- 8      # residues at that terminus used to superpose onto the resident pose
clash_cut  <- 2.5    # heavy-atom contact below this = clash, reported not fixed
seed       <- 42

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

## ---- helpers ---------------------------------------------------------------
## Kabsch superposition: rotate/translate P onto Q (paired rows, n x 3).
kabsch <- function(P, Q) {
  cP <- colMeans(P); cQ <- colMeans(Q)
  s  <- svd(t(sweep(P, 2, cP)) %*% sweep(Q, 2, cQ))
  d  <- sign(det(s$v %*% t(s$u)))          # guard against a reflection
  list(R = s$v %*% diag(c(1, 1, d)) %*% t(s$u), cP = cP, cQ = cQ)
}
apply_fit <- function(X, f) sweep(sweep(X, 2, f$cP) %*% t(f$R), 2, f$cQ, "+")

ca_xyz <- function(pdb, ch) {
  a <- pdb$atom[pdb$atom$chain == ch & pdb$atom$elety == "CA", ]
  a <- a[order(a$resno), ]
  list(resno = a$resno, xyz = as.matrix(a[, c("x", "y", "z")]))
}

min_cross <- function(A, B) {   # closest approach between two atom sets
  d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  sqrt(max(0, min(d2)))
}

AA3 <- c(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",GLU="E",GLY="G",HIS="H",
         ILE="I",LEU="L",LYS="K",MET="M",PHE="F",PRO="P",SER="S",THR="T",TRP="W",
         TYR="Y",VAL="V",SEC="U",PYL="O")

## Cleave a residue range straight out of the PDB TEXT.
##
## Deliberately does not use cleave_peptide_from_af(): that path selects on
## `elety`, so a PDB whose atom-name column holds element symbols (which a bad
## CIF->PDB conversion produces) yields zero CA atoms and reports an EMPTY
## extracted sequence -- a mismatch error that says nothing about the real
## problem. Reading columns 13-16 ourselves makes the file describe itself, and
## removes any dependence on which build of the package is installed.
##
## Residue numbering is left as-is (UniProt positions), not renumbered to 1..N.
cleave_pdb_text <- function(src, start, end, out, expected_seq = NULL) {
  raw <- grep("^(ATOM|HETATM)", readLines(src, warn = FALSE), value = TRUE)
  if (!length(raw)) stop("no atom records in ", src)

  aname <- trimws(substr(raw, 13, 16))
  rname <- trimws(substr(raw, 18, 20))
  rno   <- suppressWarnings(as.integer(substr(raw, 23, 26)))

  if (!any(aname == "CA"))
    stop("no CA atoms in ", src, "\n  atom names present: ",
         paste(head(unique(aname), 12), collapse = " "),
         "\n  resno range: ", paste(range(rno, na.rm = TRUE), collapse = "-"),
         "\n  -> this PDB has element symbols where atom names belong; delete it ",
         "so af_structure_to_pdb() rebuilds it")

  keep <- !is.na(rno) & rno >= start & rno <= end
  if (!any(keep))
    stop("no residues in ", start, "-", end, "; file covers ",
         paste(range(rno, na.rm = TRUE), collapse = "-"))

  ca  <- keep & aname == "CA"
  got <- paste0(ifelse(is.na(AA3[rname[ca]]), "X", AA3[rname[ca]])[order(rno[ca])],
                collapse = "")
  if (!is.null(expected_seq) && !identical(got, expected_seq))
    stop("sequence mismatch over ", start, "-", end,
         "\n  got:      ", got, " (", nchar(got), " res)",
         "\n  expected: ", expected_seq, " (", nchar(expected_seq), " res)")

  writeLines(c(raw[keep], "END"), out)
  list(path = out, seq = got, n_res = sum(ca))
}

## ---- 1. pick a random open conformation for this receptor ------------------
man <- readr::read_csv(manifest_csv, show_col_types = FALSE) %>%
  filter(receptor == rec_name | receptor_id == rec_id) %>%
  mutate(path = file.path(open_dir, pdb_file)) %>%
  filter(file.exists(path))

if (nrow(man) == 0)
  stop("no extracted open model for ", rec_name, " (", rec_id, ") in ", open_dir,
       "\n  -> fall back to position_ligand_initial_guess() on the apo AF-DB model")

pick <- man[sample(nrow(man), 1), ]
message(sprintf("open conformation: %s  [%s, conformer %s, opened by %s, ipTM %.3f, %d res in pocket]",
                basename(pick$path), pick$afpd_dir_name, pick$conformer,
                pick$peptide, pick$iptm, pick$num_res_in_pocket))

cx <- bio3d::read.pdb(pick$path, verbose = FALSE)

## receptor = the chain with the most residues; resident peptide = next largest
n_res  <- tapply(cx$atom$resno, cx$atom$chain, function(x) length(unique(x)))
ord    <- names(sort(n_res, decreasing = TRUE))
if (length(ord) < 2) stop("expected >=2 chains in ", pick$path, ", found ", length(ord))
rec_ch <- ord[1]; lig_ch <- ord[2]
message(sprintf("  chains: receptor %s (%d res), resident peptide %s (%d res)",
                rec_ch, n_res[[rec_ch]], lig_ch, n_res[[lig_ch]]))

## ---- 2. test peptide, cleaved fresh from AF-DB -----------------------------
lig_src <- af_structure_to_pdb(lig_id, cache_dir = af_db)
cl <- cleave_pdb_text(lig_src, lig_start, lig_end,
                      out          = file.path(out_dir, paste0(lig_label, ".pdb")),
                      expected_seq = lig_seq)
lig_pdb_path <- cl$path
lig <- bio3d::read.pdb(lig_pdb_path, verbose = FALSE)
message(sprintf("test peptide: %s  (%d res, %s)",
                basename(lig_pdb_path), cl$n_res, cl$seq))

## ---- 3. superpose the test peptide onto the resident pose ------------------
## Anchor on the pocket-binding terminus, so the residues that actually sit
## deepest in the pocket are the ones matched. The two peptides differ in
## length, so we take the last (C) or first (N) n_anchor CAs of each.
res_ca <- ca_xyz(cx,  lig_ch)
new_ca <- ca_xyz(lig, unique(lig$atom$chain)[1])
## same failure mode as the ligand file, but for the open-model complex: if its
## atom names are element symbols there are no CAs to anchor on, and Kabsch
## would fail on an empty matrix rather than say why.
if (nrow(res_ca$xyz) < 3 || nrow(new_ca$xyz) < 3)
  stop("too few CA atoms to superpose (resident ", nrow(res_ca$xyz),
       ", test ", nrow(new_ca$xyz), ")\n  elety in ", basename(pick$path), ": ",
       paste(head(unique(cx$atom$elety), 12), collapse = " "))

k <- min(n_anchor, length(res_ca$resno), length(new_ca$resno))
idx <- function(n, k) if (terminus == "C") seq.int(n - k + 1, n) else seq_len(k)

fit  <- kabsch(new_ca$xyz[idx(nrow(new_ca$xyz), k), , drop = FALSE],
               res_ca$xyz[idx(nrow(res_ca$xyz), k), , drop = FALSE])

lig_xyz_old <- as.matrix(lig$atom[, c("x", "y", "z")])
lig_xyz_new <- apply_fit(lig_xyz_old, fit)

anchor_rmsd <- sqrt(mean(rowSums((
  apply_fit(new_ca$xyz[idx(nrow(new_ca$xyz), k), , drop = FALSE], fit) -
  res_ca$xyz[idx(nrow(res_ca$xyz), k), , drop = FALSE])^2)))
message(sprintf("superposed on %d %s-terminal CAs, anchor RMSD %.2f A", k, terminus, anchor_rmsd))

lig$atom[, c("x", "y", "z")] <- lig_xyz_new
lig$xyz <- bio3d::as.xyz(as.vector(t(lig_xyz_new)))

## ---- 4. receptor with the resident peptide removed -------------------------
rec <- bio3d::trim.pdb(cx, bio3d::atom.select(cx, chain = rec_ch, verbose = FALSE))
rec_xyz <- as.matrix(rec$atom[, c("x", "y", "z")])

## ---- 5. QC: does the swapped peptide actually fit? -------------------------
closest <- min_cross(lig_xyz_new, rec_xyz)
n_clash <- sum(apply(lig_xyz_new, 1, function(p)
                 min(sqrt(colSums((t(rec_xyz) - p)^2))) < clash_cut))
message(sprintf("closest ligand-receptor approach %.2f A; %d ligand atom(s) under %.1f A",
                closest, n_clash, clash_cut))
if (n_clash > 0)
  message("  !! clashes present -- the test peptide is bulkier than the resident one here. ",
          "Try another conformer (re-run with a different seed) or let AF relax it.")

## ---- 6. write the combined initial guess -----------------------------------
## The filename stem MUST equal AlphaPulldown's job description or
## --initial_guess_dir silently skips the job (the lookup is a dict .get()).
## AP builds that description by joining one component per `+`-separated input
## part with "_and_", where a part is "<accession>" or "<accession>_<start>-<end>".
## So `--input NK2R+TKN1:72-107` -> "P21452_and_P20366_72-107".
## The conformer is recorded in the provenance CSV instead of the filename.
out_pdb <- file.path(out_dir, sprintf("%s_and_%s_%s-%s.pdb",
                                      rec_id, lig_id, lig_start, lig_end))
bio3d::write.pdb(bio3d::cat.pdb(rec, lig, rechain = TRUE), file = out_pdb)

message("initial guess: ", out_pdb)

## provenance for this build -- which open conformer it came from matters if the
## scores later need auditing for chemotype bias
readr::write_csv(
  pick %>% transmute(receptor, receptor_id, conformer, source_model = basename(path),
                     opened_by_peptide = peptide, iptm, num_res_in_pocket,
                     test_ligand = lig_label, test_ligand_id = lig_id,
                     test_range = paste0(lig_start, "-", lig_end),
                     terminus, n_anchor = k, anchor_rmsd, closest, n_clash,
                     output = basename(out_pdb)),
  sub("\\.pdb$", "_provenance.csv", out_pdb))
