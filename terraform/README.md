# Hosting

Static hosting for the per-gene LigandFinder plot pages.

    S3 (private)  <--OAC--  CloudFront  <--https--  reader

The bucket is never public. CloudFront reaches it through an Origin Access
Control, and the bucket policy admits exactly one principal: this distribution,
matched by ARN. Public access is additionally blocked at the bucket level, so
exposure would take two independent mistakes rather than one.

## Why not S3 static website hosting

Website endpoints are HTTP-only and need a public bucket. These URLs get shared
in papers and emails, so HTTPS matters; CloudFront supplies it on the default
certificate with no domain to buy.

## Why pages are pre-gzipped

Measured on the real corpus: **~15x** (a 32.5 MB page ships as 2.1 MB). Uploading
pre-compressed with `Content-Encoding: gzip` applies that saving to storage as
well as transfer, and keeps egress inside the CloudFront free tier at any
plausible traffic level. `deploy.sh` handles it.

## Interim mode: serving straight from S3

CloudFront is gated on AWS account verification for new accounts, which can take
days. `-var use_cloudfront=false` serves from the S3 REST endpoint instead:

    terraform apply -var use_cloudfront=false

Still HTTPS, and object ACLs stay blocked, but the bucket policy becomes
public-read (`s3:GetObject` only -- the bucket cannot be listed and writes still
need IAM credentials). No CDN, no custom domain, and egress comes off S3's
100 GB/month free tier rather than CloudFront's 1 TB.

Revert once the account clears -- this re-privatises the bucket and puts the
distribution in front:

    terraform apply    # use_cloudfront defaults to true

## Usage

    terraform init
    terraform apply
    ./deploy.sh ~/AF2_analysis/new_meth_plot3

If the bucket already exists (created by hand), adopt it rather than recreating:

    terraform import aws_s3_bucket.plots ligandfinder-plots-kb

## Cost

Storage is the only fixed cost: ~4.7 GB gzipped for the full ~5,200-page corpus,
about $0.11/month. Egress sits inside CloudFront's free tier unless traffic is
very high. A monthly AWS budget alarm is worth setting regardless.
