# Odds ratio
odds.ratio <- function(a, b, total){
  a.b <- length(intersect(a, b))
  a <- length(a) 
  b <- length(b) 
  a.nonb <- a - a.b
  nona.b <- b - a.b
  nona.nonb <- total - a.b - a.nonb - nona.b
  x <- matrix(c(a.b, a.nonb, nona.b, nona.nonb), 2, byrow = TRUE)
  or <- OddsRatio(x)
  or
}