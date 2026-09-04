#!/bin/bash
## ---------------------------------------------------------------------------
## Publish the per-gene ligandFinder plot pages to S3 + CloudFront.
##
## The pages are self-contained HTML with the SVG inlined, which makes them
## big (~15 MB median, 214 MB worst case, ~75 GB for the full 5177-gene set)
## but ~15x gzip-compressible.  So every object is uploaded PRE-COMPRESSED with
## Content-Encoding: gzip -- the browser inflates it transparently, and what
## crosses the wire (and what S3 bills for) is ~5 GB rather than ~75 GB.
## These are the same headers the first four test objects already carry.
##
## Files are uploaded in parallel, and every success is appended to a log, so
## a re-run only picks up what is still missing.  Interrupt it freely.
##
## Usage:
##   bash inst/scripts/upload_plots_s3.sh [src_dir] [bucket] [distribution_id]
##
##   src_dir     Directory of <gene>.html pages
##               (default: ~/AF2_analysis/ligandFinder_v3)
##   bucket      S3 bucket           (default: ligandfinder-plots-kb)
##   dist_id     CloudFront dist id  (default: E1WXP3FZJDNP71; "" = skip
##               the cache invalidation at the end)
##
## Environment:
##   JOBS=10     parallel uploads.  This is network-bound, not CPU-bound:
##               gzip -6 runs at ~110 MB/s/core, so the whole 75 GB costs
##               only ~11 core-minutes, while each job spends most of its
##               life waiting on the network.  Past ~10 you are just adding
##               aws-cli processes (~120 MB each) without filling the pipe
##               any faster.
##   GZIP_LEVEL=6
##   DRYRUN=1    compress and resolve every key, but do not PUT anything
##
## Progress prints on one live line: pages done, GB read, GB sent, current
## upload rate and ETA.  It is derived from the upload log, so it stays
## accurate across a re-run and you can also watch it from another shell:
##   wc -l < ~/AF2_analysis/ligandFinder_v3/.s3_uploaded.log
##
## Example:
##   caffeinate -i bash inst/scripts/upload_plots_s3.sh
##
## Wall time is bound by upstream bandwidth: ~5 GB of compressed payload, so
## roughly 7 min on 100 Mbit, under a minute on gigabit.  Compression itself
## runs in parallel with the uploads and is not the bottleneck.
## ---------------------------------------------------------------------------
set -uo pipefail

SRC="${1:-$HOME/AF2_analysis/ligandFinder_v3}"
BUCKET="${2:-ligandfinder-plots-kb}"
DIST="${3-E1WXP3FZJDNP71}"

JOBS="${JOBS:-10}"
GZIP_LEVEL="${GZIP_LEVEL:-6}"
DRYRUN="${DRYRUN:-0}"

SRC="${SRC%/}"
LOG="$SRC/.s3_uploaded.log"

command -v aws >/dev/null || { echo "aws cli not found" >&2; exit 1; }
[ -d "$SRC" ] || { echo "no such directory: $SRC" >&2; exit 1; }
aws sts get-caller-identity >/dev/null 2>&1 || {
  echo "aws credentials not working -- try 'aws sts get-caller-identity'" >&2; exit 1; }

TMP="$(mktemp -d "${TMPDIR:-/tmp}/lf_s3_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

## One file: compress to scratch, PUT, record the key.  Kept as a function so
## the same code path serves the dry run.
upload_one() {
  local f="$1" key rc raw gzb
  key="$(basename "$f")"
  local gz="$TMP/$key.gz"

  raw=$(stat -f %z "$f")
  gzip -"$GZIP_LEVEL" -c "$f" > "$gz" || { echo "gzip failed: $key" >&2; return 1; }
  gzb=$(stat -f %z "$gz")

  local args=(--content-type "text/html; charset=utf-8"
              --content-encoding gzip
              --cache-control "public, max-age=3600"
              --only-show-errors)
  [ "$DRYRUN" = "1" ] && args+=(--dryrun)

  aws s3 cp "$gz" "s3://$BUCKET/$key" "${args[@]}"
  rc=$?
  rm -f "$gz"

  if [ $rc -eq 0 ]; then
    ## key <TAB> raw bytes <TAB> gzipped bytes -- the byte columns are what
    ## the progress reporter integrates into a rate and an ETA.  Short lines
    ## appended O_APPEND are atomic, so parallel writers do not interleave.
    [ "$DRYRUN" = "1" ] || printf '%s\t%s\t%s\n' "$key" "$raw" "$gzb" >> "$LOG"
  else
    echo "FAILED: $key" >&2
  fi
  return $rc
}
export -f upload_one
export TMP BUCKET LOG GZIP_LEVEL DRYRUN

## Worklist = every page not already recorded as uploaded.
total=$(find "$SRC" -maxdepth 1 -name '*.html' | wc -l | tr -d ' ')
todo="$TMP/todo"
if [ -s "$LOG" ]; then
  find "$SRC" -maxdepth 1 -name '*.html' | sort > "$TMP/all"
  cut -f1 "$LOG" | sed "s|^|$SRC/|" | sort -u > "$TMP/done"
  comm -23 "$TMP/all" "$TMP/done" > "$todo"
else
  find "$SRC" -maxdepth 1 -name '*.html' | sort > "$todo"
fi

n=$(wc -l < "$todo" | tr -d ' ')
echo "source      : $SRC"
echo "bucket      : s3://$BUCKET"
echo "pages       : $n to upload ($total total, $((total - n)) already logged)"
echo "parallelism : $JOBS"
[ "$DRYRUN" = "1" ] && echo "MODE        : dry run, nothing will be written"
[ "$n" -eq 0 ] && { echo "nothing to do"; exit 0; }

## Total raw bytes across the WHOLE set, so a resumed run still reports
## progress against the real finish line rather than against what is left.
total_raw=$(find "$SRC" -maxdepth 1 -name '*.html' -exec stat -f %z {} + \
            | awk '{s+=$1} END{print s+0}')

## Where this run started, so the rate reflects this run and not an average
## dragged down by a session that was interrupted hours ago.
base_raw=0
[ -s "$LOG" ] && base_raw=$(awk -F'\t' '{s+=$2} END{print s+0}' "$LOG")

human_time() {                       # seconds -> 45s / 12m30s / 1h04m
  local t=$1
  if   [ "$t" -lt 60 ];   then printf '%ds' "$t"
  elif [ "$t" -lt 3600 ]; then printf '%dm%02ds' $((t/60)) $((t%60))
  else                         printf '%dh%02dm' $((t/3600)) $((t%3600/60)); fi
}

report() {
  local n=0 r=0 g=0 el rate eta
  if [ -s "$LOG" ]; then
    set -- $(awk -F'\t' '{n++; r+=$2; g+=$3} END{printf "%d %.0f %.0f", n+0, r+0, g+0}' "$LOG")
    n=$1; r=$2; g=$3
  fi
  el=$(( $(date +%s) - start )); [ "$el" -lt 1 ] && el=1
  rate=$(( (r - base_raw) / el ))                       # raw bytes/s this run
  if [ "$rate" -gt 0 ] && [ "$total_raw" -gt "$r" ]; then
    eta=$(human_time $(( (total_raw - r) / rate )))
  else
    eta="--"
  fi
  printf '%s  %d/%d pages  %.1f/%.1f GB read  %.2f GB sent  %.0f MB/s  ETA %s' \
    "$1" "$n" "$total" \
    "$(echo "$r"         | awk '{print $1/1073741824}')" \
    "$(echo "$total_raw" | awk '{print $1/1073741824}')" \
    "$(echo "$g"         | awk '{print $1/1073741824}')" \
    "$(echo "$rate"      | awk '{print $1/1048576}')" "$eta"
}

start=$(date +%s)

## Live line while the uploads run.  On a terminal it rewrites one line; piped
## to a file (nohup, tee) it prints a fresh line every 30s so the log stays
## readable instead of filling with carriage returns.
PROG_PID=""
if [ "$DRYRUN" != "1" ]; then
  if [ -t 1 ]; then
    while :; do report $'\r'; sleep 5; done &
  else
    while :; do report ''; echo; sleep 30; done &
  fi
  PROG_PID=$!
  trap 'rm -rf "$TMP"; [ -n "$PROG_PID" ] && kill "$PROG_PID" 2>/dev/null' EXIT
fi

tr '\n' '\0' < "$todo" | xargs -0 -P "$JOBS" -I{} bash -c 'upload_one "$@"' _ {}
rc=$?

if [ -n "$PROG_PID" ]; then kill "$PROG_PID" 2>/dev/null; wait "$PROG_PID" 2>/dev/null; fi
[ "$DRYRUN" != "1" ] && { report ''; echo; }
echo "elapsed: $(human_time $(( $(date +%s) - start )))"

if [ "$DRYRUN" != "1" ] && [ -s "$LOG" ]; then
  left=$(( total - $(cut -f1 "$LOG" | sort -u | wc -l | tr -d ' ') ))
  if [ "$left" -gt 0 ]; then
    echo "$left page(s) still missing -- re-run this script to retry just those" >&2
  fi
fi

## The pages are served with max-age=3600, so anything already in the edge
## cache (the four test objects) would otherwise serve stale for up to an hour.
## A single "/*" path counts as one invalidation, not one per object.
if [ "$DRYRUN" != "1" ] && [ $rc -eq 0 ] && [ -n "$DIST" ]; then
  echo "invalidating CloudFront cache ($DIST) ..."
  aws cloudfront create-invalidation --distribution-id "$DIST" --paths '/*' \
    --query 'Invalidation.{Id:Id,Status:Status}' --output text
fi

exit $rc
