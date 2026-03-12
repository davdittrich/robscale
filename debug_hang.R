#!/usr/bin/env Rscript
library(robscale)
set.seed(42)

cat("Testing qn(n=128)...\n")
x <- rnorm(128)
q <- qn(x)
cat("qn(n=128) OK, result:", q, "\n")

cat("Testing sn(n=128)...\n")
s <- sn(x)
cat("sn(n=128) OK, result:", s, "\n")

cat("Testing robLoc(n=128)...\n")
l <- robLoc(x)
cat("robLoc(n=128) OK, result:", l, "\n")

cat("Testing qn(n=16384)...\n")
x2 <- rnorm(16384)
q2 <- qn(x2)
cat("qn(n=16384) OK, result:", q2, "\n")

cat("Testing sn(n=16384)...\n")
s2 <- sn(x2)
cat("sn(n=16384) OK, result:", s2, "\n")

cat("Testing robLoc(n=16384)...\n")
l2 <- robLoc(x2)
cat("robLoc(n=16384) OK, result:", l2, "\n")

cat("All isolated tests passed.\n")
