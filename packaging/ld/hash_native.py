import hashlib
from pathlib import Path
import sys
p = Path(sys.argv[1])
p.with_suffix(p.suffix + '.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest() + '\n')
