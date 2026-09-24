# Validation environment notes

Windows R 4.6.1 uses the private landgraph 0.0.3 library. RStudio's bundled
Pandoc is used for all vignette builds. Locale variables inherited as C.UTF-8
are cleared before starting R.

At approximately 17:31 on September 24, the concurrent source-test and covr
jobs each had 24 threads on a 20-logical-processor machine, at 100% CPU load.
The first source-test run was interrupted before it returned test results.
Its only emitted line was `Testing terradish`; it is not passing evidence.
Moving that log failed because the shell still held it open. The replacement
run then reused the log path, so there is no separate original log file.

The replacement source-test run and final check set OMP_NUM_THREADS,
OPENBLAS_NUM_THREADS, and MKL_NUM_THREADS to 1. The vignette job finished cleanly.
Coverage had left instrumented .o/.dll/.gc* artifacts in the shared candidate,
so the replacement source tests loaded those artifacts. That source-test
attempt was also interrupted before results, with its log retained as
phase8_tests_instrumented_interrupted.log. Ordinary tests now use the separate
phase8-source-tests copy, built from clean sources.

The initial coverage process still had its large thread pool, with substantial
CPU use and no completed coverage receipt after about thirty minutes. It was
stopped and restarted in its own phase8-coverage copy with the same one-thread
limits. The initial runner and log are retained with initial/interrupted names.
The final archive check was already using clean source files in its own build
directory. Only completed replacement receipts establish final validation.
These environment controls do not change the package API or stopping criteria.

The restarted default native-coverage run reached gcov collection but failed
with exit code 6 while processing amg_solver.gcno. It produced no usable
coverage data, and is not counted as passing coverage. The error is retained
in phase8_coverage.log/rds. A separate covr::package_coverage run measures all
R-source lines under the full test suite with covr.flags=character() and
covr.gcov="", preserving normal native compilation. This matches the plan's
Phases 2-5 target, which consists of R files. Its results must be labeled R
coverage; C++ numerical tests pass independently, but native line coverage is
not established by the failed collector.

The R-only coverage run completed successfully: 86.79266% overall and
88.02621% weighted coverage across the 25 R files changed in Phases 2-5.
Three individual changed files remain below 80%; RELEASE_REVIEW.md lists them.
Final source tests, installed tests, and the final archive check also completed.
No interrupted run is used as evidence for these results.

The final check runs through Rscript with vignette rebuilding enabled. A
literal interactive RStudio check has not been performed. The same RStudio
Pandoc executable is used, but this does not establish a GUI-session check.
