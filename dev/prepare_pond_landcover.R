library(terra)
raw <- rast('dev/amphibian_example/inputs_context_30600m/nlcd_2023_30m.tif')
data_env <- new.env()
load('data/ohio_ponds.rda', envir = data_env)
template <- unwrap(data_env$ohio_ponds$covariates)[[1]]
# A dominant native class is nominal. It differs from mean forest fraction.
nominal <- resample(raw, template, method = 'mode')
cover <- nominal
values(cover) <- ifelse(values(nominal)[, 1] %in% c(41, 42, 43, 90), 1L, 2L)
cover <- as.factor(cover)
levels(cover) <- data.frame(ID = 1:2, class = c('forest', 'nonforest'))
names(cover) <- 'landcover'
stopifnot(!anyNA(values(cover)), compareGeom(template, cover, stopOnError = FALSE))
saveRDS(wrap(cover, proxy = FALSE), 'inst/extdata/ohio-ponds/landcover100.rds', version = 2)
print(freq(cover))
