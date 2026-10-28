"""Stage package files only; preserve the shared checkout and its libraries."""
from pathlib import Path
import shutil

root = Path(__file__).resolve().parent.parent
validation = root / "dev/individual_vignette_integration_20261023"
stage = validation / "package"
stage.mkdir(parents=True, exist_ok=True)
for name in ("R", "data", "inst", "man", "src", "tests", "vignettes"):
    shutil.copytree(root / name, stage / name, dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("*_cache", "*_files", "*.o", "*.dll"))
for name in ("DESCRIPTION", "NAMESPACE", "LICENSE", "NEWS.md", ".Rbuildignore",
             "README.md", "configure", "configure.ac", "configure.in",
             "configure.win", "cleanup", "cleanup.win"):
    if (root / name).is_file():
        shutil.copy2(root / name, stage / name)
(validation / "library").mkdir(exist_ok=True)
print(stage)
