files <- list.files("man", pattern="[.]Rd$",full.names=TRUE)
issues <- lapply(files, function(f) tools::checkRd(tools::parse_Rd(f)))
names(issues) <- basename(files)
bad <- issues[lengths(issues)>0L]
print(bad)
stopifnot(length(bad)==0L)
cat(length(files), "Rd files passed checkRd.\n")
