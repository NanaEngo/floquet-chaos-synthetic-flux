#!/usr/bin/env python3
"""Hash reconstruction source files for a reproducibility manifest."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

def digest(path: Path) -> str:
    h=hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest()

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument('--root', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a=p.parse_args(); root=a.root.resolve(); files=sorted(root.glob('scripts/*.py'))
    data={'root':str(root), 'files':{str(f.relative_to(root)):digest(f) for f in files}}
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8'); print(json.dumps(data,indent=2,sort_keys=True)); return 0
if __name__=='__main__': raise SystemExit(main())
