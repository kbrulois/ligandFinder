# The MMseqs2 search settings for the candidate pipeline, in ONE place: sourced
# by candidate_pipeline.sh (which runs them) and candidate_manifest.sh (which
# records them), so the record cannot drift from what ran.
#
# Defaults find nothing on peptides this short: the k-mer prefilter is tuned
# for proteins, and a candidate has to share an exact k-mer with a target
# before it is aligned at all -- so 8-10 aa candidates got NO hit, and only a
# fifth of 11-15 aa ones did (2026-10-27). --prefilter-mode 1 replaces the k-mer
# stage with an ungapped alignment of the candidate against EVERY target, which
# the 8-60 aa target set (~670 sequences) makes cheap; --max-seqs 100 then
# passes the best 100 of those to the gapped alignment. Checked against no
# prefilter at all (--prefilter-mode 2) on 2,000 candidates: the best hit was
# identical for 100%, the best hit within each target type for 99.95%.
# (--exhaustive-search 1 is NOT the same thing: it found fewer hits than the
# k-mer prefilter, not more.) The e-value is opened right up so every
# candidate keeps its closest known peptide, however weak; judge a hit by its
# e-value / bits_per_res, not by its presence. Composition bias correction is
# off: it is tuned for proteins and penalises short peptides.
MMSEQS_ARGS=(--prefilter-mode 1 --max-seqs 100 --min-seq-id 0 -e 10000 --comp-bias-corr 0
             --alignment-mode 3 --format-output "query,target,pident,alnlen,bits,evalue,qcov,tcov")
