#!/usr/bin/env python3
"""Fail closed on prohibited PZHEEV references; PZHEEVD is allowed."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys

SYMBOL = re.compile(r'(?<![A-Za-z0-9_])_?pzheev_?(?![A-Za-z0-9_])', re.I)
SUFFIXES = {'.f', '.f90', '.f95', '.f03', '.for', '.c', '.cc', '.cpp', '.h', '.hpp'}

def fortran_code(text, fixed=False, literals=False):
    result = []
    for line in text.splitlines():
        if fixed and line and line[0] in 'cC*!':
            continue
        out = []; quote = None; i = 0
        while i < len(line):
            c = line[i]
            if quote:
                if literals: out.append(c)
                if c == quote:
                    if i + 1 < len(line) and line[i+1] == quote:
                        if literals: out.append(line[i+1])
                        i += 2; continue
                    quote = None
            elif c in "\"'":
                quote = c
                if literals: out.append(c)
            elif c == '!':
                break
            else:
                out.append(c)
            i += 1
        result.append(''.join(out))
    return re.sub(r'&[ \t]*\n[ \t]*&?', '', '\n'.join(result))

def strict_files(root):
    def fail(error):
        raise error
    for directory, dirs, files in os.walk(str(root), onerror=fail):
        for name in dirs:
            if Path(directory, name).is_symlink():
                raise RuntimeError('directory symlink cannot be verified: ' + str(Path(directory, name)))
        for name in sorted(files):
            yield Path(directory, name)

def dynamic_reference(code):
    # Tokenize strings as indivisible tokens so documentation cannot look like code.
    tokens = re.findall(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_][A-Za-z_0-9]*|[()]', code)
    for i, token in enumerate(tokens[:-1]):
        if token.lower() not in {'bind', 'dlsym'} or tokens[i+1] != '(':
            continue
        depth = 0
        for item in tokens[i+1:]:
            if item == '(': depth += 1
            elif item == ')':
                depth -= 1
                if depth == 0: break
            elif item[:1] in {'"', "'"} and SYMBOL.fullmatch(item[1:-1]):
                return True
    return False

def source_findings(root):
    root = Path(root)
    if not root.is_dir():
        raise RuntimeError('source root does not exist: ' + str(root))
    hits = []
    for p in strict_files(root):
        if p.suffix.lower() not in SUFFIXES or not p.is_file():
            continue
        if getattr(p.stat(), 'st_flags', 0) & 0x40000000:
            raise RuntimeError('cloud-only source cannot be verified: ' + str(p))
        text = p.read_text(encoding='utf-8')
        if p.suffix.lower() in {'.f', '.f90', '.f95', '.f03', '.for'}:
            fixed = p.suffix.lower() in {'.f', '.for'}
            preserved = fortran_code(text, fixed, literals=True)
            code = fortran_code(text, fixed)
        else:
            lexical = r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\''
            preserved = re.sub(lexical, lambda m: '' if m.group().startswith(('/',)) else m.group(), text, flags=re.S)
            code = re.sub(lexical, '', text, flags=re.S)
        bindings = dynamic_reference(preserved)
        if SYMBOL.search(code) or bindings:
            hits.append(str(p))
    return hits

def object_findings(paths):
    hits = []
    for p in paths:
        result = subprocess.run(['nm', '-u', str(p)], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode:
            raise RuntimeError('cannot inspect object ' + str(p) + ': ' + result.stderr.strip())
        if any(SYMBOL.fullmatch(line.split()[-1]) for line in result.stdout.splitlines() if line.split()):
            hits.append(str(p))
    return hits

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--object-tree', type=Path)
    parser.add_argument('--built-executable', type=Path, help='new build output to make non-executable if inspection fails')
    args = parser.parse_args()
    if not args.source_root and not args.object_tree:
        parser.error('provide --source-root or --object-tree')
    def quarantine():
        if args.built_executable and args.built_executable.exists():
            args.built_executable.chmod(args.built_executable.stat().st_mode & ~0o111)
    try:
        hits = source_findings(args.source_root) if args.source_root else []
        if args.object_tree:
            if not args.object_tree.is_dir():
                raise RuntimeError('object tree does not exist: ' + str(args.object_tree))
            objects = sorted(p for p in strict_files(args.object_tree) if p.suffix.lower() in {'.o', '.obj'})
            if not objects:
                raise RuntimeError('no application objects found: ' + str(args.object_tree))
            hits += object_findings(objects)
        if hits:
            quarantine()
            print('ERROR: PZHEEV is prohibited. Use a validated PZHEEVD implementation.\n' + '\n'.join(hits), file=sys.stderr)
            return 1
        print('PZHEEV prohibition check PASS')
        return 0
    except (OSError, UnicodeError, RuntimeError) as exc:
        quarantine()
        print('ERROR: prohibition check could not verify safety: ' + str(exc), file=sys.stderr)
        return 2

if __name__ == '__main__':
    sys.exit(main())
