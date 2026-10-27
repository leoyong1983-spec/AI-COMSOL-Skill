"""Check distributable integrity and run solver-independent behavioral tests."""
import ast
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "skills/ai-comsol-skill"


def main():
    required = [ROOT/"README.md", ROOT/"LICENSE", ROOT/"THIRD_PARTY_NOTICES.md",
                SKILL/"SKILL.md", SKILL/"agents/openai.yaml", ROOT/"docs/verification.md"]
    for path in required:
        if not path.is_file():
            raise ValueError("missing required file: " + str(path.relative_to(ROOT)))
    count=0
    for parent in (SKILL, ROOT/"tools", ROOT/"docs"):
        for path in parent.rglob("*"):
            if not path.is_file() or path.suffix not in (".md", ".py", ".cs", ".json", ".yaml"):
                continue
            text=path.read_bytes().decode("utf-8")
            if "\ufffd" in text:
                raise ValueError("replacement character: " + path.name)
            if path.suffix==".json":
                json.loads(text, parse_constant=lambda x: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
            if path.suffix==".py":ast.parse(text, filename=str(path))
            if path.suffix==".md":
                for target in re.findall(r"\]\(([^)]+)\)",text):
                    if "://" in target or target.startswith("#") or target.startswith("mailto:"):continue
                    dest=target.split("#",1)[0]
                    if dest and not (path.parent/dest).exists():
                        raise ValueError("broken local link in " + path.name + ": " + target)
            count+=1
    entry=(SKILL/"SKILL.md").read_text(encoding="utf-8")
    if not entry.startswith("---\nname: ai-comsol-skill\n"):
        raise ValueError("skill name/frontmatter mismatch")
    if "$ai-comsol-skill" not in (SKILL/"agents/openai.yaml").read_text(encoding="utf-8"):
        raise ValueError("UI prompt does not invoke the skill")
    for name in ("parameter-sources.json", "catalyst-reactor-data-v2.json"):
        data=json.loads((SKILL/"references"/name).read_text(encoding="utf-8"))
        for source in data["sources"]:
            if not all(source.get(k) for k in ("id","grade","access","license","relation_to_comsol_runtime","adoption_decision")):
                raise ValueError("incomplete source provenance in " + name)
    if "MIT License" not in (ROOT/"LICENSE").read_text(encoding="utf-8"):
        raise ValueError("MIT license missing")
    print(json.dumps({"status":"repository-integrity-passed","text_files_checked":count}))
    return subprocess.call([sys.executable,"-B","-m","unittest","discover","-s",str(SKILL/"scripts"),"-p","test_*.py","-v"],cwd=ROOT)


if __name__=="__main__":
    raise SystemExit(main())
