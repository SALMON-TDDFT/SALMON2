"""Synchronize embedded guide text with Markdown; --check is read-only."""
import argparse
import json
from pathlib import Path
import re


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]/'salmonpy'
    def read(name):
        base = root.parent if name.startswith('docs/') else root
        return (base/name).read_text(encoding='utf-8')
    path=root/'manual.json'
    original=path.read_text(encoding='utf-8')
    data=json.loads(original)
    for name in data['documents']:
        text=read(name)
        data['documents'][name]={'sections':[{'heading':name,'text':text}]}
    for task in data['tasks']:
        guides=task.get('guides') or {'guide':task['guide']}
        task['instructions']=[{'heading':language,'text':read(name)}
                              for language,name in guides.items()]
        if task['id'] in ('public_api','snapshots'):
            task['examples']=list(dict.fromkeys(code for name in guides.values()
                for code in re.findall(r'```python\n(.*?)```',read(name),re.S)))
    result=json.dumps(data,ensure_ascii=False,indent=2)+'\n'
    if args.check:
        if result!=original: raise SystemExit('manual.json guide text is out of sync; run tools/sync_manual_guides.py')
    else:path.write_text(result,encoding='utf-8')


if __name__=='__main__':main()
