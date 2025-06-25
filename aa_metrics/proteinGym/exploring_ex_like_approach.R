rm(list = ls(all.names = TRUE))

a = 0.278
b = 0.666

cumul_power <- function(t) { (t ** b) / (a + t ** b) }

inverse_cumul_power <- function(p) {
  return((p * a / (1 - p)) ^ (1/b))
}


cumul_dist2 <- sapply(1:100 / 100, cumul_power)

inverse_cumul_power(0.47)

pdf_power <- function(t) {
  numerator <- b * a * (t^(b-1))
  denominator <- (a + t^b)^2
  return(numerator / denominator)
}

pdf_power(0.12)