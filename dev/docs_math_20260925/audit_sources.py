"""Check active documentation encodings and inventory Pandoc math nodes."""
from pathlib import Path
import json
import re
import subprocess

root = Path(__file__).resolve().parents[2]
out = Path(__file__).resolve().parent
pandoc = Path("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe")
files = [Path(n) for n in subprocess.check_output(
    ["git", "ls-files"], cwd=root, text=True).splitlines()]
files = [p for p in files if p.parts[0] not in {"dev", "src"}
         and p.suffix in {".md", ".Rmd", ".Rd", ".R", ".txt"}]
issues = []
for rel in files:
    for i, line in enumerate((root / rel).read_text(encoding="utf-8-sig").splitlines(), 1):
        if re.search("[\u00c3\u00c2\u00e2\u00cf\u00ce\ufffd]|[\x00-\x08\x0b\x0c\x0e-\x1f]", line):
            issues.append({"file": str(rel), "line": i, "text": line})

def math_nodes(node):
    if isinstance(node, dict):
        if node.get("t") == "Math":
            yield {"type": node["c"][0]["t"], "tex": node["c"][1]}
        for value in node.values():
            yield from math_nodes(value)
    elif isinstance(node, list):
        for value in node:
            yield from math_nodes(value)

inventory = {}
for p in sorted((root / "vignettes").glob("*.Rmd")):
    text = p.read_text(encoding="utf-8")
    assert not re.search(r"[\\/]d?eqn\{", text), p
    # R Markdown chunk headers are interpreted by knitr, not raw Pandoc.
    # Remove fenced code so R's $ operator cannot become a false math node.
    text = re.sub(r"(?ms)^```[^\n]*\n.*?^```[^\n]*(?:\n|$)", "\n", text)
    parsed = subprocess.run([str(pandoc), "--from=markdown", "--to=json"],
                            input=text, text=True, encoding="utf-8",
                            capture_output=True, check=True)
    inventory[p.name] = list(math_nodes(json.loads(parsed.stdout)))
record = {"files_scanned": len(files), "encoding_issues": issues,
          "vignette_math": inventory}
(out / "source_inventory.json").write_text(
    json.dumps(record, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
assert not issues, issues
print(f"Encoding scan: {len(files)} files, no suspicious sequences.")
for name, equations in inventory.items():
    print(f"{name}: {len(equations)} parsed math expressions")
