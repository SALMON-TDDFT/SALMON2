"""Integration assertion for the restored native HSE Lowdin path."""
import re,sys,math
from pathlib import Path
s=Path(sys.argv[1]).read_text()
a=re.findall(r'HSE_PT_CN_ORTH pre_gram=\s*(\S+) post_gram=\s*(\S+) pre_residual=\s*(\S+) post_residual=\s*(\S+) correction=\s*(\S+)',s)
assert a, 'Missing native HSE Lowdin hook diagnostics'
for row in a:
 p,q,rb,ra,c=map(float,row)
 assert all(math.isfinite(v) for v in [p,q,rb,ra,c])
 assert p<=1e-8 and q<1e-12 and rb<1e-10
assert 'end SALMON' in s
print('LOWDIN_INTEGRATION_PASS',len(a))
