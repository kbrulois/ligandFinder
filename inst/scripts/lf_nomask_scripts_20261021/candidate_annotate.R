## Stage 3 of the candidate pipeline: fold the MMseqs2 search back onto the
## candidate table and add the derived columns (Kevin, 2026-10-23).
##
##   HEAD=xgb_pm THRESH=0.5 TOL=4 Rscript candidate_annotate.R
##
## Hits are named <gene>x<start>x<end>, the same shape as a candidate
## peptide_id, so a hit reads directly against a candidate. The readable name
## and the accession are SEPARATE columns: "Activation peptide" is the UniProt
## description on 77 different accessions, so a name is not an identifier and
## must not be grouped on.
##
## Columns added
##   mmseqs_best_hit       <gene>x<start>x<end> of the closest known peptide
##   mmseqs_best_hit_name  its readable name   mmseqs_best_hit_acc  its accession
##   mmseqs_hit_type       ligand_list (a curated bioactive ligand) | peptide
##                         (a UniProt peptide feature) | propeptide (a real
##                         cleavage product but NOT a signalling ligand) | none
##   mmseqs_pident/bits/evalue/qcov/tcov   of that hit. Best is chosen by BIT
##                         SCORE -- by pident a 6-residue exact match would beat
##                         a 25-residue match at 90%. tcov is needed to read
##                         pident: 100% over a 20-residue slice of a 1,400-residue
##                         chain is not a rediscovery; pident 100 AND tcov ~1 is.
##   mmseqs_bits_per_res   bits / candidate length; raw bits grow with length and
##                         cannot be compared across sizes
##   ll_* / pep_* / prop_* the best hit WITHIN each class. 40% of the target set
##                         is propeptides, so the overall best hit alone conflates
##                         "looks like a hormone" with "looks like a pro-segment";
##                         a candidate topped by a propeptide still shows its best
##                         curated-ligand match here.
##   length_bin            5-10, 11-20, ... to LEN_MAX
##   overlap_*             whether the candidate's RANGE intersects an annotation
##                         on its OWN precursor -- a different question from the
##                         mmseqs_* columns, which find a similar sequence
##                         anywhere. A candidate can resemble a peptide on another
##                         gene and overlap nothing here, or overlap a known
##                         propeptide and resemble nothing.
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
HEAD    <- Sys.getenv("HEAD", "xgb_pm")
thr     <- as.numeric(Sys.getenv("THRESH", "0.5"))
TOL     <- as.integer(Sys.getenv("TOL", "4"))
LEN_MAX <- as.integer(Sys.getenv("LEN_MAX", "60"))
O       <- file.path(A2, "mmseq_cand")
f_cand  <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_T%s_tol%d.csv", HEAD, thr, TOL))

cd  <- fread(f_cand)
kn  <- fread(file.path(O, "known.csv"))
ids <- fread(file.path(O, "cand_ids.csv"))
m   <- fread(file.path(O, "res.m8"), header = FALSE,
             col.names = c("qid","target","pident","alnlen","bits","evalue","qcov","tcov"))

## Carry the target's identity and the query's sequence onto every hit with
## UPDATE joins (X[Y, on=, col := i.col]): they add columns and can never
## multiply rows, unlike X[Y, on=], which fans out on any duplicated key.
m[kn, on = c(target = "tid"), `:=`(ttype = i.type, kid = i.known_id,
                                   kname = i.name, kacc = i.accession)]
m[ids, on = "qid", sequence := i.sequence]
m <- m[!is.na(ttype) & !is.na(sequence)]

pfx  <- c(ligand_list = "ll", peptide = "pep", propeptide = "prop")
## wipe any previous annotation so re-running is idempotent
cols <- c(paste0("mmseqs_", c("best_hit","best_hit_name","best_hit_acc","hit_type",
                              "pident","bits","evalue","qcov","tcov","bits_per_res")),
          as.vector(outer(unname(pfx), c("_hit","_hit_name","_pident","_tcov","_bits_per_res"),
                          paste0)),
          "length_bin", "overlap_known","overlap_id","overlap_name","overlap_type",
          "overlap_frac_cand","overlap_exact")
if (length(intersect(names(cd), cols))) cd[, (intersect(names(cd), cols)) := NULL]

## ---- best hit overall ------------------------------------------------------
setorder(m, sequence, -bits, evalue, kid)
best <- m[m[, .I[1], by = sequence]$V1]
cd[best, on = "sequence", `:=`(
     mmseqs_best_hit = i.kid, mmseqs_best_hit_name = i.kname, mmseqs_best_hit_acc = i.kacc,
     mmseqs_hit_type = i.ttype, mmseqs_pident = i.pident, mmseqs_bits = i.bits,
     mmseqs_evalue = i.evalue, mmseqs_qcov = i.qcov, mmseqs_tcov = i.tcov)]
cd[is.na(mmseqs_best_hit), `:=`(mmseqs_best_hit = "none", mmseqs_best_hit_name = "none",
                                mmseqs_best_hit_acc = "none", mmseqs_hit_type = "none")]
cd[, mmseqs_bits_per_res := round(mmseqs_bits / length, 4)]

## ---- best hit within each class -------------------------------------------
setorder(m, sequence, ttype, -bits, evalue, kid)
per <- m[m[, .I[1], by = .(sequence, ttype)]$V1]
for (tt in names(pfx)) {
  p <- pfx[[tt]]
  cd[per[ttype == tt], on = "sequence",
     (paste0(p, c("_hit","_hit_name","_pident","_tcov","_bits_per_res"))) :=
       .(i.kid, i.kname, i.pident, i.tcov, round(i.bits / length, 4))]
  cd[is.na(get(paste0(p, "_hit"))), (paste0(p, c("_hit","_hit_name"))) := .("none", "none")]
}

## ---- length bins -----------------------------------------------------------
brk <- c(4, 10, seq(20, ceiling(LEN_MAX / 10) * 10, by = 10))
brk <- unique(brk[brk <= max(c(LEN_MAX, cd$length))])
if (max(cd$length) > max(brk)) brk <- c(brk, max(cd$length))
cd[, length_bin := cut(length, breaks = brk, labels = paste0(head(brk, -1) + 1, "-", brk[-1]),
                       include.lowest = FALSE)]
stopifnot(!any(is.na(cd$length_bin)))

## ---- positional overlap with annotations on the SAME precursor -------------
kr <- kn[!is.na(start) & !is.na(end) & end >= start,
         .(accession, k_start = start, k_end = end, k_id = known_id, k_name = name, k_type = type)]
q <- cd[, .(i = .I, accession, c_start = start, c_end = end, c_len = length)]
j <- kr[q, on = "accession", allow.cartesian = TRUE, nomatch = 0]
j[, shared := pmin(k_end, c_end) - pmax(k_start, c_start) + 1L]
j <- j[shared > 0]
setorder(j, i, -shared, k_id)
bv <- j[j[, .I[1], by = i]$V1]
cd[, `:=`(overlap_known = FALSE, overlap_id = "none", overlap_name = "none",
          overlap_type = "none", overlap_frac_cand = 0, overlap_exact = FALSE)]
cd[bv$i, `:=`(overlap_known = TRUE, overlap_id = bv$k_id, overlap_name = bv$k_name,
              overlap_type = bv$k_type, overlap_frac_cand = round(bv$shared / bv$c_len, 4),
              overlap_exact = bv$c_start == bv$k_start & bv$c_end == bv$k_end)]

fwrite(cd, f_cand)
options(width = 230)
cat(sprintf("annotated %s rows, %d columns -> %s\n", format(nrow(cd), big.mark = ","),
            ncol(cd), basename(f_cand)))
cat("\nbest hit by class:\n")
print(cd[, .(candidates = .N, pct = round(100 * .N / nrow(cd), 1)), by = mmseqs_hit_type][order(-candidates)],
      row.names = FALSE)
cat("\nbest hit WITHIN each class:\n")
print(data.table(class = names(pfx),
                 with_hit = sapply(pfx, function(p) sum(cd[[paste0(p, "_hit")]] != "none")),
                 median_pident = sapply(pfx, function(p) round(median(cd[[paste0(p, "_pident")]], na.rm = TRUE), 1))),
      row.names = FALSE)
cat(sprintf("\noverlapping an annotation on their own precursor: %s (%.0f%%), exact range match: %s\n",
            format(sum(cd$overlap_known), big.mark = ","), 100 * mean(cd$overlap_known),
            format(sum(cd$overlap_exact), big.mark = ",")))
print(cd[overlap_known == TRUE, .N, by = overlap_type][order(-N)], row.names = FALSE)
