"""AI-facing offline task search: python -m salmonpy.manual バンド"""
import json
from importlib.resources import files


def load():
    return json.loads(files('salmonpy').joinpath('manual.json').read_text(encoding='utf-8'))


def search(query):
    """Return matching task records, with executable guidance and source metadata."""
    if not isinstance(query, str) or not query.strip():
        raise ValueError('query must be nonempty text')
    data = load()
    tokens = query.casefold().split()
    ranked = []
    for task in data['tasks']:
        phrases = [task['id'], task['intent'], *task['keywords'], *task['api']]
        score = sum(1 for token in tokens for phrase in phrases
                    if token in phrase.casefold() or phrase.casefold() in token)
        if score:
            ranked.append((score, task))
    return {'schema_version':data['schema_version'], 'query':query,
            'tasks':[t for _, t in sorted(ranked,key=lambda x:-x[0])],
            'getting_started':data['getting_started'],
            'task_catalog':[{'id':t['id'], 'intent':t['intent']} for t in data['tasks']],
            'global_rules':data['global_rules'], 'sources':data['sources'],
            'environment':data['environment'], 'beginner_intake':data['beginner_intake']}


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('query',nargs='+')
    args=parser.parse_args()
    print(json.dumps(search(' '.join(args.query)),ensure_ascii=False,indent=2))
