## Which MSA each precursor gets, and from where (Kevin, 2026-11-10: "replace
## the existing conservation channel with conservation scores computed from
## alphafold db msa").
##
## AFDB's own MSAs have answered anonymous callers with 403 since Aug 2026, so
## the alignments are ColabFold/MMseqs2 ones (UniRef30 2023-02 + the ColabFold
## env DB 2021-08): the LiteFold human-proteome set on Hugging Face wherever it
## holds a precursor's sequence EXACTLY, the ColabFold MMseqs2 API for the rest.
## The match is by SEQUENCE, never by accession: 16 of the 5,243 cache
## accessions the set shares carry another isoform or sequence revision there,
## and an alignment whose query row is not the cache sequence describes other
## residues (msa_conservation.py refuses one).
##
## This writes into D the manifest and the three lists the fetchers and the
## scorer read:
##   manifest.tsv         the LiteFold manifest (acc, ..., sequence, a3m_path,
##                        ..., a3m_bytes), downloaded unless MANIFEST exists
##   hf_needed.tsv        acc, a3m_bytes -- ONE file per distinct cache sequence
##                        the set holds (its first accession in manifest order),
##                        for fetch_hf_msas.sh. So 7 precursors have no file of
##                        their own: AMY1A and AMY1C read AMY1B's, OPN1MW2 and
##                        OPN1MW3 read OPN1MW's, IFNA1 reads P0DY56's, ...
##   api_needed.tsv       accession, seq -- the precursors whose sequence the set
##                        does not hold, for fetch_colabfold_api.py
##   hf_needed_alias.tsv  acc, sequence of every hf_needed file, for
##                        msa_conservation.py --alias. The scorer finds those 7
##                        files by sequence, among the cache's accessions and
##                        these; P0DY56 (the set's IFNA1, P01562 in the cache)
##                        is the one only the alias names
##
## 2026-11-10 (manifest 12,139,597 bytes, md5 ae209eb51389eca51f38f07782e7648d,
## 18,519 entries): 5,616 precursors, 5,602 distinct sequences; 5,228 held by
## the set -> 5,222 files, 30.9 GB of A3M; 388 for the API (303,169 residues,
## the longest 14,507 -- immunoglobulin segments, sequence revisions, ...).
## The manifest lives on `main`, so a fresh download can differ; point MANIFEST
## at a saved copy to rebuild those lists exactly.
##
## Then, from inst/scripts/generate_residue_db/msa_cons/:
##   D=$D P=2 ./fetch_hf_msas.sh                 # HF rate-limits; P=8 drew 429s
##   ~/.virtualenvs/lf-plm/bin/python -I fetch_colabfold_api.py $D/api_needed.tsv $D/api_colabfold --workers 3
##   ~/.virtualenvs/lf-plm/bin/python msa_conservation.py --workers 4 [--allow-missing]
## msa_conservation.py's --msa-dirs, --alias, --out and --per-protein default to
## the default D; for any other, pass --msa-dirs $D/hf_colabfold $D/api_colabfold
## --alias $D/hf_needed_alias.tsv and an --out / --per-protein of your choosing.
##
##   LC_ALL=en_US.UTF-8 Rscript inst/scripts/generate_residue_db/msa_cons/msa_fetch_lists.R
##   D=/tmp/msa MANIFEST=~/AF2_analysis/lf_msa/manifest.tsv LC_ALL=en_US.UTF-8 Rscript msa_fetch_lists.R
##
## D        output directory (default ~/AF2_analysis/lf_msa, fetch_hf_msas.sh's D)
## MANIFEST the manifest to read (default $D/manifest.tsv); downloaded from
##          MANIFEST_URL when absent, and copied into D when it is elsewhere,
##          because fetch_hf_msas.sh reads $D/manifest.tsv
## CACHE    the precursor list (default ~/AF2_analysis/lf_pepend_residue_cache.rds,
##          the peptide-end residue cache; precursor_seqs.tsv is its $prec)
suppressMessages(library(data.table))
D        <- path.expand(Sys.getenv("D", "~/AF2_analysis/lf_msa"))
MANIFEST <- path.expand(Sys.getenv("MANIFEST", file.path(D, "manifest.tsv")))
URL      <- Sys.getenv("MANIFEST_URL",
  "https://huggingface.co/datasets/LiteFold/human-proteome-wide-msa/resolve/main/manifest.tsv")
CACHE    <- path.expand(Sys.getenv("CACHE", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
dir.create(D, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(MANIFEST)) {
  ## to a temporary name first: a cut transfer must not pass for a manifest
  tmp <- paste0(MANIFEST, ".tmp")
  rc <- system2("curl", c("-sfL", "--retry", "3", "-o", shQuote(tmp), shQuote(URL)))
  if (!(rc == 0 && file.exists(tmp) && isTRUE(startsWith(readLines(tmp, n = 1L), "acc\t")))) {
    unlink(tmp); stop("manifest download failed (curl rc ", rc, "): ", URL, call. = FALSE)
  }
  file.rename(tmp, MANIFEST)
  message("downloaded ", URL)
}
if (normalizePath(MANIFEST) != normalizePath(file.path(D, "manifest.tsv"), mustWork = FALSE))
  invisible(file.copy(MANIFEST, file.path(D, "manifest.tsv"), overwrite = TRUE))
cat(sprintf("manifest %s: %s bytes, md5 %s\n", MANIFEST, format(file.size(MANIFEST), big.mark = ","),
            unname(tools::md5sum(MANIFEST))))

p <- readRDS(CACHE)$prec
cat("precursors", nrow(p), "| distinct accessions", length(unique(p$accession)),
    "| distinct sequences", length(unique(p$seq)), "\n")
m <- as.data.frame(fread(MANIFEST, select = c("acc", "sequence", "a3m_bytes")))
cat("manifest entries", nrow(m), "\n")

hit_seq <- p$seq %in% m$sequence
cat("cache sequences held by the set, exactly:", sum(hit_seq), "\n")
j <- merge(data.frame(accession = p$accession, seq = p$seq), m, by.x = "accession", by.y = "acc")
cat("accessions shared with the set", nrow(j), "| of those, same sequence", sum(j$seq == j$sequence), "\n")

## one file per distinct sequence: the first manifest row that carries it
mm <- m[m$sequence %in% p$seq, ]
mm <- mm[!duplicated(mm$sequence), ]
cat("HF files", nrow(mm), "| GB", sum(mm$a3m_bytes) / 1e9, "\n")
need <- p[!hit_seq, ]
cat("for the ColabFold API:", nrow(need), "| residues", sum(nchar(need$seq)), "| longest", max(nchar(need$seq)), "\n")
print(summary(nchar(need$seq)))

fwrite(mm[, c("acc", "a3m_bytes")], file.path(D, "hf_needed.tsv"), sep = "\t")
fwrite(data.frame(accession = need$accession, seq = need$seq), file.path(D, "api_needed.tsv"), sep = "\t")
## every manifest row of an hf_needed accession, in manifest order
fwrite(m[m$acc %in% mm$acc, c("acc", "sequence")], file.path(D, "hf_needed_alias.tsv"), sep = "\t")
cat("wrote", file.path(D, c("hf_needed.tsv", "api_needed.tsv", "hf_needed_alias.tsv")), sep = "\n  ")
