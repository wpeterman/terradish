# Curate the validated, frozen example for the installed package.
# Run from the package root. This does not fit models or regenerate genotypes.
bundle <- file.path("dev", "amphibian_example",
  "terrestrial_individual_candidate_20261023", "vignette-candidate")
stopifnot(file.exists(file.path(bundle, "DATASET_READINESS.json")))
ready <- jsonlite::fromJSON(file.path(bundle, "DATASET_READINESS.json"))
stopifnot(isTRUE(ready$ready_for_focused_SNP_IBR_vignette))
destination <- file.path("inst", "extdata", "individual-salamanders")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

individual_salamanders <- readRDS(file.path(bundle, "data", "individual_salamanders.rds"))
save(individual_salamanders, file = "data/individual_salamanders.rda",
     compress = "xz", version = 2)

# Keep every original prediction, including unsuccessful marker comparisons.
paths <- list.files(file.path(bundle, "results"), "__predictions[.]csv$", full.names = TRUE)
predictions <- setNames(lapply(paths, read.csv),
  sub("__predictions[.]csv$", "", basename(paths)))
stopifnot(length(predictions) == 40L)
saveRDS(predictions, file.path(destination, "predictions.rds"), compress = "xz", version = 2)
assessment <- jsonlite::fromJSON(file.path(bundle, "results", "SUMMARY.json"))
assay <- jsonlite::fromJSON(file.path(bundle, "results", "ASSAY_SUMMARY.json"))
validation <- jsonlite::fromJSON(file.path(bundle, "results", "VALIDATION_SUMMARY.json"))
results <- list(original = assessment$comparisons,
  independent_assay = validation$comparisons,
  block_contrasts = assay$contrast_central_95_ranges,
  independent_scores = read.csv(file.path(bundle, "results", "VALIDATION_SCORES.csv")))
saveRDS(results, file.path(destination, "assessment.rds"), compress = "xz", version = 2)

# The second assay remains available for checking the frozen models.
independent <- lapply(names(individual_salamanders$data), function(arm) {
  read_panel <- function(suffix) {
    tab <- read.csv(file.path(bundle, "data", paste0(arm, suffix)), check.names = FALSE)
    ans <- as.matrix(tab[-1L])
    rownames(ans) <- tab$sample_id
    ans
  }
  list(snp = read_panel("_validation_snp_dosage.csv"),
       microsat = read_panel("_validation_microsat_alleles.csv"))
})
names(independent) <- names(individual_salamanders$data)
saveRDS(independent, file.path(destination, "independent-assay.rds"), compress = "xz", version = 2)
for (name in c("drying_terms_10m.tif", "environment_10m.tif", "terrain_30m.tif"))
  stopifnot(file.copy(file.path(bundle, "data", name), file.path(destination, name), overwrite = TRUE))

# Attribution and protocol receipts are data, not a required JSON dependency.
provenance <- c("landscape_download_manifest.json", "SAMPLING_PROTOCOL.md",
  "validation_PLAN.json", "validation_MODEL_FREEZE_FOR_VALIDATION.json",
  "validation_SEED_MANIFEST.json", "ANALYSIS_SOURCE_SHA256.json")
dir.create(file.path(destination, "provenance"), showWarnings = FALSE)
for (name in provenance)
  stopifnot(file.copy(file.path(bundle, "provenance", name),
    file.path(destination, "provenance", name), overwrite = TRUE))
for (name in c("BOOTSTRAP_ESTIMATES.csv", "panel_500.json", "panel_1000.json")) {
  source <- file.path(bundle, "results", name)
  if (file.exists(source)) stopifnot(file.copy(source,
    file.path(destination, "provenance", name), overwrite = TRUE))
}
for (name in c("landscape-and-samples", "drying-response-and-markers", "predictive-comparisons"))
  stopifnot(file.copy(file.path(bundle, "figures", paste0(name, ".png")),
    file.path("vignettes", "precomputed", paste0("individual-salamanders-", name, ".png")),
    overwrite = TRUE))
message("Curated two matched scenarios, 40 prediction records, and the frozen second assay.")
