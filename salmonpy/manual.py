"""AI-facing offline task search: python -m salmonpy.manual バンド"""
import json
import re
from importlib.resources import files


def load():
    return json.loads(files('salmonpy').joinpath('manual.json').read_text(encoding='utf-8'))


def search(query, *, lang='ja'):
    """Return matching task records, with executable guidance and source metadata."""
    if not isinstance(query, str) or not query.strip():
        raise ValueError('query must be nonempty text')
    if lang not in ('en', 'ja'): raise ValueError('lang must be en or ja')
    data = load()
    stopwords = {'i', 'a', 'an', 'the', 'to', 'do', 'can', 'how', 'what', 'want',
                 'run', 'use', 'please', 'with', 'of', 'in', 'is', 'my'}
    tokens = [t for t in re.findall(r'[a-z0-9_]+|[^\x00-\x7f]+', query.casefold())
              if t not in stopwords]
    def matches(token, phrase):
        if token.isascii():
            return re.search(r'(?<![a-z0-9_])'+re.escape(token)+r'(?![a-z0-9_])', phrase.casefold()) is not None
        return token in phrase.casefold() or phrase.casefold() in token
    ranked = []
    for task in data['tasks']:
        phrases = [task['id'], task['intent'], *task['keywords'], *task['api']]
        score = sum(1 for token in tokens for phrase in phrases if matches(token, phrase))
        if task['id'].casefold() in tokens: score += 100
        if score:
            localized = dict(task)
            localized['intent'] = task.get('intent_i18n',{}).get(lang,task['intent'])
            ranked.append((score, localized))
    return {'schema_version':data['schema_version'], 'query':query,
            'tasks':[t for _, t in sorted(ranked,key=lambda x:-x[0])],
            'language':lang,
            'getting_started':data.get('getting_started_i18n',{}).get(lang,data['getting_started']),
            'task_catalog':[{'id':t['id'], 'intent':t.get('intent_i18n',{}).get(lang,t['intent'])} for t in data['tasks']],
            'global_rules':data['global_rules'], 'sources':data['sources'],
            'environment':data['environment'], 'beginner_intake':data.get('beginner_intake_i18n',{}).get(lang,data['beginner_intake'])}


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lang', choices=('en','ja'), default='ja')
    parser.add_argument('query',nargs='+')
    args=parser.parse_args()
    print(json.dumps(search(' '.join(args.query), lang=args.lang),ensure_ascii=False,indent=2))
