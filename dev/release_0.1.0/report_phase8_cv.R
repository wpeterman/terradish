x <- readRDS("dev/release_0.1.0/phase8_cv_modes_comparison.rds")
x$selected <- x$extended_minus_base > 1e-6
x$tied <- abs(x$extended_minus_base) <= 1e-6
print(aggregate(x[c("selected", "tied", "failures", "baseline_failures")],
                x[c("family", "nuisance")], sum))
fixed <- x[x$nuisance == "fixed", ]
stopifnot(all(tapply(fixed$selected, fixed$family, sum) <= 3))
write.csv(x, "dev/release_0.1.0/phase8_cv_modes_comparison.csv", row.names = FALSE)
