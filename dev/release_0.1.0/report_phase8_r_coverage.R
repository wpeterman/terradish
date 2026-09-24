receipt <- readRDS("dev/release_0.1.0/phase8_r_coverage.rds")
coverage <- receipt$result
stopifnot(inherits(coverage, "coverage"))
rows <- covr::tally_coverage(coverage, by = "line")
rows$filename <- gsub("\\", "/", rows$filename, fixed = TRUE)
files <- lapply(split(rows, rows$filename), function(x)
  data.frame(file = x$filename[1], lines = nrow(x), covered = sum(x$value > 0),
             percent = 100 * mean(x$value > 0)))
files <- do.call(rbind, files)
rownames(files) <- NULL
touched <- system2("git", c("--no-optional-locks", "diff", "604a3a1", "9b75ff4",
  "--name-only", "--diff-filter=ACMR", "--", "R"), stdout = TRUE)
files$touched_phases2to5 <- files$file %in% touched
write.csv(files, "dev/release_0.1.0/phase8_r_coverage_by_file.csv", row.names = FALSE)
write.csv(rows[rows$value == 0 & rows$filename %in% touched, ],
  "dev/release_0.1.0/phase8_r_uncovered_touched_lines.csv", row.names = FALSE)
print(files)
target <- files[files$touched_phases2to5, ]
stopifnot(nrow(target) == length(touched))
cat("Overall line coverage:", covr::percent_coverage(coverage), "\n")
cat("Changed R-file weighted line coverage:",
    100 * sum(target$covered) / sum(target$lines), "\n")
cat("Changed files below 80%:\n")
print(target[target$percent < 80, ])
