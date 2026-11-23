"""AI-facing offline task search: python -m salmonpy.manual バンド"""
import json
import re
from importlib.resources import files


def load():
    return json.loads(files('salmonpy').joinpath('manual.json').read_text(encoding='utf-8'))


def search(query, *, lang='ja', limit=5):
    """Return bounded task summaries; retrieve one guide with task_details."""
    if not isinstance(query, str) or not query.strip():
        raise ValueError('query must be nonempty text')
    if lang not in ('en', 'ja'): raise ValueError('lang must be en or ja')
    if type(limit) is not int or limit < 1:
        raise ValueError('limit must be a positive integer')
    data = load()
    stopwords = {'i', 'a', 'an', 'the', 'to', 'do', 'can', 'how', 'what', 'want',
                 'run', 'use', 'please', 'with', 'of', 'in', 'is', 'my'}
    tokens = [t for t in re.findall(r'[a-z0-9_]+|[^\x00-\x7f]+', query.casefold())
              if t not in stopwords]
    def matches(token, phrase):
        token=token.casefold()
        if token.isascii():
            return re.search(r'(?<![a-z0-9_])'+re.escape(token)+r'(?![a-z0-9_])', phrase.casefold()) is not None
        return token in phrase.casefold() or phrase.casefold() in token
    ranked = []
    for task in data['tasks']:
        phrases = [task['id'], task['intent'], *task.get('intent_i18n',{}).values(),
                   *task['keywords'], *task['api']]
        score = sum(1 for token in tokens for phrase in phrases if matches(token, phrase))
        if task['id'].casefold() in tokens: score += 100
        components=task.get('workflow_components',())
        hits=[component for component in components if matches(component,query)]
        required=task.get('workflow_required',())
        # A bare list such as "TDDFT FDTD" is also a combined-workflow query.
        # Additional words (e.g. "compare") must not imply coupling by default.
        remaining=query.casefold()
        for component in hits:
            remaining=re.sub(r'(?<![a-z0-9_])'+re.escape(component.casefold())+
                             r'(?![a-z0-9_])','',remaining)
        models_only=not re.sub(r'\band\b|と|[\s,+/&]+','',remaining)
        if (len(hits)>=2 and any(component in hits for component in required)
                and (models_only or any(word.casefold() in query.casefold()
                                           for word in task.get('workflow_keywords',())))):
            score += 100 * len(hits)
        if score:
            localized = {k:task[k] for k in ('id','api','inputs','outputs','constraints','environment',
                                             'guide','guides','sample_files') if k in task}
            localized['intent'] = task.get('intent_i18n',{}).get(lang,task['intent'])
            ranked.append((score, localized))
    return {'schema_version':data['schema_version'], 'query':query,
            'tasks':[t for _, t in sorted(ranked,key=lambda x:-x[0])[:limit]],
            'matches':len(ranked),
            'details_command':'python -m salmonpy.manual --lang '+lang+' --task TASK_ID',
            'language':lang,
            'getting_started':data.get('getting_started_i18n',{}).get(lang,data['getting_started']),
            'task_catalog':[{'id':t['id'], 'intent':t.get('intent_i18n',{}).get(lang,t['intent'])} for t in data['tasks']],
            'global_rules':data['global_rules'], 'sources':data['sources'],
            'environment':data['environment'], 'beginner_intake':data.get('beginner_intake_i18n',{}).get(lang,data['beginner_intake'])}


def task_details(task_id, *, lang='ja'):
    """Return one exact task with only the requested guide language.

    A task without a translated guide returns its declared source guide, with
    instruction_language identifying that fallback rather than relabeling it.
    """
    if lang not in ('en','ja'):
        raise ValueError('lang must be en or ja')
    task=next((t for t in load()['tasks'] if t['id']==task_id),None)
    if task is None:
        raise ValueError('unknown task ID: '+str(task_id))
    result=dict(task)
    result['intent']=task.get('intent_i18n',{}).get(lang,task['intent'])
    instructions=task['instructions']
    selected=next((item for item in instructions if item['heading']==lang),instructions[0])
    result['instructions']=[selected]
    result['instruction_language']=selected['heading'] if selected['heading'] in ('ja','en') else 'source'
    result['translation_available']=selected['heading']==lang
    result['requested_language']=lang
    return result


def main():
    import sys
    import argparse
    if hasattr(sys.stdout,'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lang',choices=('en','ja'),default='ja')
    parser.add_argument('--task',help='exact task ID; return this task and its guide only')
    parser.add_argument('--limit',type=int,default=5,help='maximum search summaries (default5)')
    parser.add_argument('query',nargs='*')
    args=parser.parse_args()
    if args.task and args.query:
        parser.error('choose either a query or --task TASK_ID')
    try:
        if args.task:
            result=task_details(args.task,lang=args.lang)
        elif args.query:
            result=search(' '.join(args.query),lang=args.lang,limit=args.limit)
        else:
            result=load()['getting_started_i18n'][args.lang]
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
