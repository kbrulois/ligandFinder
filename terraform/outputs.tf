output "bucket" {
  description = "Upload target for deploy.sh"
  value       = aws_s3_bucket.plots.id
}

output "distribution_id" {
  description = "Cache invalidation target; null while serving straight from S3"
  value       = one(aws_cloudfront_distribution.plots[*].id)
}

output "site_url" {
  description = "Base URL; per-gene pages are <site_url>/<GENE>.html"
  value       = var.use_cloudfront ? "https://${one(aws_cloudfront_distribution.plots[*].domain_name)}" : "https://${aws_s3_bucket.plots.id}.s3.${var.region}.amazonaws.com"
}

output "example_deep_link" {
  description = "Deep link format the pages resolve via their data-peps attributes"
  value       = var.use_cloudfront ? "https://${one(aws_cloudfront_distribution.plots[*].domain_name)}/ANO8.html#ANO8_w702-737" : "https://${aws_s3_bucket.plots.id}.s3.${var.region}.amazonaws.com/ANO8.html#ANO8_w702-737"
}
