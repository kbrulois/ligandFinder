variable "region" {
  description = "AWS region for the bucket. CloudFront itself is global."
  type        = string
  default     = "us-west-2"
}

variable "bucket_name" {
  description = "S3 bucket holding the generated per-gene pages. Globally unique across all of AWS."
  type        = string
  default     = "ligandfinder-plots-kb"
}

variable "price_class" {
  # PriceClass_100 = North America + Europe edge locations only. The full class
  # (PriceClass_All) adds Asia/South America edges at higher per-GB rates; not
  # worth it for a site whose readers are overwhelmingly US/EU academics.
  description = "CloudFront edge coverage."
  type        = string
  default     = "PriceClass_100"
}

variable "use_cloudfront" {
  # CloudFront is gated on AWS account verification for new accounts, which can
  # take days. false serves straight from the S3 REST endpoint instead -- still
  # HTTPS, but the bucket must be publicly readable and there is no CDN or
  # custom domain. Intended as a documented interim; flip back to true once the
  # account clears and the distribution takes over with the bucket private again.
  description = "Serve via CloudFront (private bucket) rather than a public S3 bucket."
  type        = bool
  default     = true
}
