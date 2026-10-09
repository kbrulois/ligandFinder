# LigandFinder plot hosting: private S3 bucket fronted by CloudFront.
#
# Why not an S3 static-website bucket? Website endpoints are HTTP-only and
# require the bucket itself to be public. CloudFront gives HTTPS on the default
# *.cloudfront.net certificate and lets the bucket stay fully private -- the
# distribution reaches it through an Origin Access Control, and nothing else can.
#
# State is local (terraform.tfstate on disk). For a single operator that is fine;
# a shared setup would put state in S3 with DynamoDB locking so two people cannot
# apply at once.

terraform {
  required_version = ">= 1.5"
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  region = var.region
}

resource "aws_s3_bucket" "plots" {
  bucket = var.bucket_name
}

# Belt and braces: even if a future policy or ACL tried to expose the bucket,
# these four flags override it. Access is meant to arrive only via CloudFront.
resource "aws_s3_bucket_public_access_block" "plots" {
  bucket = aws_s3_bucket.plots.id

  # ACLs stay blocked in both modes -- nothing here should ever be granted via
  # an object ACL. Only the two POLICY-related flags relax, and only when
  # serving publicly, because a public-read bucket policy is precisely what they
  # exist to prevent.
  block_public_acls       = true
  ignore_public_acls      = true
  block_public_policy     = var.use_cloudfront
  restrict_public_buckets = var.use_cloudfront
}

# OAC is the current mechanism (it signs origin requests with SigV4). It replaces
# Origin Access Identity, which is legacy and does not support all regions/KMS.
resource "aws_cloudfront_origin_access_control" "plots" {
  count = var.use_cloudfront ? 1 : 0

  name                              = "${var.bucket_name}-oac"
  description                       = "CloudFront -> private S3 origin"
  origin_access_control_origin_type = "s3"
  signing_behavior                  = "always"
  signing_protocol                  = "sigv4"
}

data "aws_cloudfront_cache_policy" "optimized" {
  name = "Managed-CachingOptimized"
}

resource "aws_cloudfront_distribution" "plots" {
  count = var.use_cloudfront ? 1 : 0

  enabled     = true
  comment     = "LigandFinder per-gene plot pages"
  price_class = var.price_class

  origin {
    domain_name              = aws_s3_bucket.plots.bucket_regional_domain_name
    origin_id                = "s3-${var.bucket_name}"
    origin_access_control_id = aws_cloudfront_origin_access_control.plots[0].id
  }

  default_cache_behavior {
    target_origin_id       = "s3-${var.bucket_name}"
    viewer_protocol_policy = "redirect-to-https"
    allowed_methods        = ["GET", "HEAD"]
    cached_methods         = ["GET", "HEAD"]
    cache_policy_id        = data.aws_cloudfront_cache_policy.optimized.id

    # Objects are uploaded ALREADY gzipped, with Content-Encoding: gzip set at
    # upload time (deploy.sh). CloudFront only compresses responses that lack
    # that header, so leaving this off simply states the intent rather than
    # relying on that behaviour. Pre-compressing also means the ~15x saving
    # applies to S3 storage, not just to transfer.
    compress = false
  }

  # Everything under the gated prefix goes through a viewer-request function
  # that demands HTTP basic auth. The rest of the bucket -- the full per-gene
  # set at the root -- keeps the open default behaviour above, so existing
  # links are unaffected. One prompt covers the whole prefix: browsers cache
  # basic credentials per origin and realm, so the gallery and every gene page
  # under it are a single login.
  dynamic "ordered_cache_behavior" {
    for_each = var.gated_prefix == "" ? [] : [1]
    content {
      path_pattern           = "${var.gated_prefix}/*"
      target_origin_id       = "s3-${var.bucket_name}"
      viewer_protocol_policy = "redirect-to-https"
      allowed_methods        = ["GET", "HEAD"]
      cached_methods         = ["GET", "HEAD"]
      cache_policy_id        = data.aws_cloudfront_cache_policy.optimized.id
      compress               = false

      function_association {
        event_type   = "viewer-request"
        function_arn = aws_cloudfront_function.basic_auth[0].arn
      }
    }
  }

  restrictions {
    geo_restriction {
      restriction_type = "none"
    }
  }

  viewer_certificate {
    cloudfront_default_certificate = true
  }
}

# The only way into the bucket: this distribution, identified by ARN. Without the
# SourceArn condition any CloudFront distribution in any AWS account could read
# it. Note this policy is not "public" in S3's sense, so it coexists with the
# public access block above.
data "aws_iam_policy_document" "cloudfront_only" {
  count = var.use_cloudfront ? 1 : 0

  statement {
    sid       = "AllowCloudFrontServicePrincipalReadOnly"
    actions   = ["s3:GetObject"]
    resources = ["${aws_s3_bucket.plots.arn}/*"]

    principals {
      type        = "Service"
      identifiers = ["cloudfront.amazonaws.com"]
    }

    condition {
      test     = "StringEquals"
      variable = "AWS:SourceArn"
      values   = [aws_cloudfront_distribution.plots[0].arn]
    }
  }
}

# Interim: anyone may GET objects. Write access is unaffected -- it still
# requires the IAM user's credentials. Scoped to GetObject only, so the bucket
# cannot be listed.
data "aws_iam_policy_document" "public_read" {
  count = var.use_cloudfront ? 0 : 1

  statement {
    sid       = "PublicReadGetObject"
    actions   = ["s3:GetObject"]
    resources = ["${aws_s3_bucket.plots.arn}/*"]

    principals {
      type        = "*"
      identifiers = ["*"]
    }
  }
}

resource "aws_s3_bucket_policy" "plots" {
  bucket = aws_s3_bucket.plots.id
  policy = var.use_cloudfront ? data.aws_iam_policy_document.cloudfront_only[0].json : data.aws_iam_policy_document.public_read[0].json
}


# ---------------------------------------------------------------------------
# HTTP basic auth at the edge, for the gated prefix only.
#
# A CloudFront Function cannot read a secret at runtime -- there is no
# environment, no KMS, no network. The credential is therefore COMPILED INTO
# the function source below and is visible to anyone who can read this AWS
# account or the terraform state. That is inherent to edge basic auth; it is a
# gate against casual access, not a secret-bearing auth system. Rotate it by
# changing var.gated_auth and re-applying.
#
# `authorization` is one of the few headers CloudFront Functions see on a
# viewer request, which is what makes this work without Lambda@Edge (and
# without its cold starts or its us-east-1 requirement).
# ---------------------------------------------------------------------------
resource "aws_cloudfront_function" "basic_auth" {
  count   = (var.use_cloudfront && var.gated_prefix != "") ? 1 : 0
  name    = "${var.bucket_name}-basic-auth"
  runtime = "cloudfront-js-2.0"
  comment = "HTTP basic auth for /${var.gated_prefix}/*"
  publish = true

  code = <<-JS
    function handler(event) {
      var want = "Basic ${base64encode(var.gated_auth)}";
      var got  = event.request.headers.authorization;
      if (got && got.value === want) return event.request;
      return {
        statusCode: 401,
        statusDescription: "Unauthorized",
        headers: {
          "www-authenticate": { value: "Basic realm=\"LigandFinder\", charset=\"UTF-8\"" },
          "cache-control":    { value: "no-store" }
        }
      };
    }
  JS
}
