from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
assert 'EXX_SR canonical occupied source; MLWF bypassed' in s,'canonical SR route not executed'
assert 'EXX_SPATIAL refresh/iterations/status/spread/gradient/overlap' not in s,'unexpected MLWF work'
assert 'EXX_SUPPORT_ACE accepted: F' not in s
import re
assert re.search(r'EXX_SR active/radius[^\n]*:  T',s),'active SR path not used'
