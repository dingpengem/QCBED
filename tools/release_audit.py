"""Read-only repository/documentation audit; writes JSON to stdout."""
from __future__ import annotations
import ast
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]


def main():
    counts = {}
    for name in ('src', 'tests', 'examples', 'reference_data', 'legacy_fortran',
                 'docs', 'benchmarks', 'tools', 'packaging'):
        paths = [p for p in (ROOT / name).rglob('*') if p.is_file()
                 and '__pycache__' not in p.parts]
        counts[name] = {'files': len(paths), 'bytes': sum(p.stat().st_size for p in paths)}
    broken = []
    docs = [ROOT / 'README.md', ROOT / 'PROJECT_SPEC.md', *sorted((ROOT / 'docs').rglob('*.md')),
            *sorted((ROOT / 'packaging').rglob('*.md'))]
    for path in docs:
        for target in re.findall(r'\[[^\]\n]+\]\(([^)\n]+)\)', path.read_text()):
            target = target.strip('<>').split('#')[0]
            if not target or re.match(r'^[a-zA-Z]+:', target):
                continue
            if not (path.parent / unquote(target)).exists():
                broken.append({'file': str(path.relative_to(ROOT)), 'target': target})
    users = ['README.md', 'docs/INSTALLATION.md', 'docs/QUICK_START.md',
             'docs/user_manual.md', 'docs/EXAMPLE_FCC_AL.md', 'examples/fcc_al.yml']
    users = sorted(set(users) | {str(p.relative_to(ROOT)) for p in (ROOT / 'examples').rglob('*')
                                if p.suffix in ('.yml', '.yaml')})
    absolute_paths = [p for p in users if re.search(r'/Users/|/home/', (ROOT / p).read_text())]
    syntax = []
    for directory in ('src', 'tests', 'tools', 'packaging'):
        for p in (ROOT / directory).rglob('*.py'):
            try:
                ast.parse(p.read_text(), filename=str(p))
            except SyntaxError as exc:
                syntax.append(str(exc))
    tracked = subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()
    generated = [p for p in tracked if any(x in Path(p).parts for x in ('__pycache__', '.pytest_cache', 'build', 'dist', 'tmp'))]
    result = dict(inventory=counts, broken_local_links=broken, user_absolute_paths=absolute_paths,
                  syntax_errors=syntax, tracked_generated_paths=generated)
    print(json.dumps(result, indent=2))
    return int(bool(broken or absolute_paths or syntax or generated))


if __name__ == '__main__':
    raise SystemExit(main())
