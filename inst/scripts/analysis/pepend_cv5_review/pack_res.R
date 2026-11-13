ROOT <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
source(file.path(ROOT, "R", "dcnn_bridge.R"))
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache_allres.rds")
prec <- fc$prec; feats <- fc$feats
stopifnot(identical(names(feats), prec$accession),
          all(vapply(feats, nrow, integer(1)) == nchar(prec$seq)))
message(sprintf("packing %d precursors, %s residues x %d channels ...",
                length(feats), format(sum(nchar(prec$seq)), big.mark=","), length(fc$all_params3)))
lf_dcnn_pack_residues(prec, feats, fc$all_params3,
                      path.expand("~/AF2_analysis/lf_pepend_residues_allres.npz"),
                      py_path = file.path(ROOT, "inst", "python"))
message("done")
