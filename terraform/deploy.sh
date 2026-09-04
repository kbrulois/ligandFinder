#!/usr/bin/env bash
# Upload generated plot pages to the CloudFront-backed bucket.
#
#   ./deploy.sh ~/AF2_analysis/new_meth_plot3
#
# Pages are gzipped before upload and stored with Content-Encoding: gzip, so the
# ~15x saving applies to S3 storage as well as to transfer. The object key drops
# the .gz, keeping URLs as <GENE>.html.
set -euo pipefail

SRC="${1:?usage: deploy.sh <directory of .html>}"
BUCKET="$(terraform -chdir="$(dirname "$0")" output -raw bucket)"
DIST="$(terraform -chdir="$(dirname "$0")" output -raw distribution_id 2>/dev/null || true)"

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

n=0
for f in "$SRC"/*.html; do
  name="$(basename "$f")"
  gzip -9 -c "$f" > "$TMP/$name.gz"
  aws s3 cp "$TMP/$name.gz" "s3://$BUCKET/$name" \
    --content-encoding gzip \
    --content-type "text/html; charset=utf-8" \
    --cache-control "public, max-age=3600" \
    --only-show-errors
  n=$((n + 1))
done
echo "uploaded $n page(s) to s3://$BUCKET/"

# max-age=3600 means edges would otherwise serve stale pages for up to an hour
# after a redeploy. Invalidating everything is free for the first 1,000 paths a
# month; beyond that, invalidate only the genes that changed.
# Skipped when serving straight from S3 (use_cloudfront = false): there is no
# edge cache, so uploads are live immediately.
if [ -n "$DIST" ] && [ "$DIST" != "null" ]; then
  aws cloudfront create-invalidation --distribution-id "$DIST" --paths '/*' \
    --query 'Invalidation.{Id:Id,Status:Status}' --output table
else
  echo "no CloudFront distribution -- objects are live immediately"
fi
