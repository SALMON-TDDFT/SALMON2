#!/usr/bin/env python3
import argparse, subprocess
from collections import defaultdict
p=argparse.ArgumentParser();p.add_argument('--binary',required=True);a=p.parse_args()
cases=[(164,4096,4096,1,27648,True),(164,4096,512,1,1280,True),
       (8,8,3,1,8,True),(8,8,3,2,16,True),(8,8,3,2,52,True),(8,8,3,1,17,True),(8,8,3,2,25,True),
       (3,10,3,1,7,True),(9,10,3,1,8,False),(1,3,3,1,4,False),
       (1,1024,1048576,2,1,False),(1,3,3,1,1,True),(0,3,3,1,1,False)]
for n,ng,nk,nc,np,ok in cases:
 out=subprocess.check_output([a.binary],input=f'{n} {ng} {nk} {nc} {np}\n',text=True)
 lines=[list(map(int,l.split())) for l in out.splitlines()];err,groups,chunk,peak,active=lines[0]
 assert (err==0)==ok,(n,ng,np,lines[0])
 if err:continue
 rows=defaultdict(list);groupdata={};work=[]
 assert len(lines[1:])==np
 for rank,r,c,count,g,leader in lines[1:]:
  assert 0<=rank<np and 1<=r<=n and count>0 and 0<=g<groups
  assert 0<=leader<np
  rows[r].append((c,c+count));work.append(count)
  if g in groupdata: assert groupdata[g]==(c,count,leader)
  groupdata[g]=(c,count,leader)
  assert lines[leader+1][0]==leader and lines[leader+1][4]==g
 assert len(rows)==n and len(groupdata)==groups
 for intervals in rows.values():
  end=1
  for lo,hi in sorted(intervals):assert lo==end;end=hi
  assert end==ng+1
 assert sum(work)==n*ng
 # Each row is split as evenly as its rank fanout permits. Different row
 # fanouts can necessarily produce a larger global spread (e.g. 3x10/7).
 fanouts=[len(intervals) for intervals in rows.values()]
 assert max(fanouts)-min(fanouts)<=1
 for intervals in rows.values():
  sizes=[hi-lo for lo,hi in intervals]
  assert max(sizes)-min(sizes)<=1
 assert active==np and peak==chunk==max(work)
 assert chunk*nk*16*nc*nc<=67108864
 print('BALANCED_COVERAGE_PASS',n,ng,np,min(work),max(work),groups)
print('BALANCED_ALL_PASS',len(cases))
