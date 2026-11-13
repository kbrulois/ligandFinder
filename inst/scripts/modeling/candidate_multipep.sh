#!/bin/zsh
# Stage 4: MultiPep class probabilities for the candidate sequences.
#
# Caches by SEQUENCE in ~/AF2_analysis/multipep_ext/multipep_scores.csv, so a
# pipeline re-run at another THRESH/TOL only scores sequences it has not seen.
# MultiPep wants ONE BARE SEQUENCE PER LINE (not FASTA -- it reads ">" headers
# as sequences), 2-200 residues, the 20 standard amino acids only.
set -u
S=${0:A:h}; A2=~/AF2_analysis; M=$A2/multipep_ext
PY=~/.virtualenvs/multipep/bin/python
REPO=$M/MultiPep
CACHE=$M/multipep_scores.csv

[[ -x $PY ]] || { echo "no multipep venv at $PY -- see ALLSEC/CANDIDATES README"; exit 1; }
[[ -d $REPO ]] || { echo "no MultiPep clone at $REPO"; exit 1; }

LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript -e '
suppressMessages(library(data.table))
A2 <- path.expand("~/AF2_analysis"); M <- file.path(A2,"multipep_ext")
u <- fread(file.path(A2,"mmseq_cand","cand_ids.csv"))
ok <- !grepl("[^ARNDCQEGHILKMFPSTWYV]", u$sequence) & nchar(u$sequence) >= 2 & nchar(u$sequence) <= 200
u <- u[ok]
f <- file.path(M,"multipep_scores.csv")
done <- if (file.exists(f)) fread(f)$sequence else character(0)
todo <- setdiff(u$sequence, done)
writeLines(todo, file.path(M,"todo_in.txt"))
cat(sprintf("  %s sequences, %s cached, %s to score\n", format(nrow(u),big.mark=","),
            format(length(done),big.mark=","), format(length(todo),big.mark=",")))' 2>&1 | grep -vE "built under|^Warning"

n=$(wc -l < $M/todo_in.txt | tr -d ' ')
if [[ $n -gt 0 ]]; then
  cd $REPO
  $PY MultiPep_predict.py -input_file $M/todo_in.txt -output_file $M/todo_out.csv -pred_type mean \
    > $M/multipep_run.log 2>&1 || { echo "MultiPep FAILED, see $M/multipep_run.log"; exit 1; }
  LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript -e '
  suppressMessages(library(data.table))
  M <- file.path(path.expand("~/AF2_analysis"),"multipep_ext")
  new <- fread(file.path(M,"todo_out.csv")); setnames(new, 1, "sequence")
  f <- file.path(M,"multipep_scores.csv")
  all <- if (file.exists(f)) unique(rbind(fread(f), new, fill=TRUE), by="sequence") else new
  fwrite(all, f); cat(sprintf("  cache now %s sequences\n", format(nrow(all),big.mark=",")))' 2>&1 | grep -vE "built under|^Warning"
else
  echo "  nothing new to score"
fi
rm -f $M/todo_in.txt $M/todo_out.csv
LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/candidate_multipep_merge.R
