#!/usr/bin/env python3
"""Fail-closed manuscript lint for unsupported scientific claims."""
from __future__ import annotations
import argparse, re
from pathlib import Path

FORBIDDEN_ACTIVE = [
    r'50\\?x', r'10--100\\?x', r'quantum advantage', r'Wigner negativity',
    r'QFI grows exponentially', r'topological protection', r'exceptional-point sensing',
]

def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument('manuscript', type=Path); p.add_argument('--require-draft-marker', action='store_true'); a=p.parse_args()
    text=a.manuscript.read_text(encoding='utf-8')
    errors=[]
    for pattern in FORBIDDEN_ACTIVE:
        for m in re.finditer(pattern, text, flags=re.I):
            context=text[max(0,m.start()-220):m.end()+220]
            negative = context.lower()
            is_withdrawn = any(marker in negative for marker in (
                'no claim', 'does not', 'do not', 'not made', 'not a ',
                'withdrawn', 'earlier claims', 'cannot be obtained', 'without direct evidence',
            ))
            if 'pending' not in negative and not is_withdrawn:
                errors.append((pattern, context.replace('\n',' ')))
    if a.require_draft_marker and 'Internal draft' not in text and 'NOT SUBMISSION READY' not in text:
        errors.append(('missing draft status',''))
    if errors:
        for pattern, context in errors: print(f'FAIL {pattern}: {context}')
        return 1
    print('PASS: no unsupported active claim detected in draft')
    return 0
if __name__=='__main__': raise SystemExit(main())
