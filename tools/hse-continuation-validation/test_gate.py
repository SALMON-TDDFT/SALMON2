from pathlib import Path
import tempfile,subprocess,json,math,sys,shutil
w=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as td:
 p=Path(td);cfg={'start_step':12536,'nt':12544,'exact_interval':16,'mpi':160,'nodes':40};(p/'preflight.json').write_text(json.dumps(cfg));(p/'launcher-exit.txt').write_text('0')
 d=p/'output.fixture/0/1';d.mkdir(parents=True)
 for rank in range(160):
  ident={'rank':rank,'omp':'12','affinity':list(range((rank%4)*12,(rank%4+1)*12)),'host':'h%d'%(rank//4)}
  (p/('summary.rank-%d.json'%rank)).write_text(json.dumps({'identity':ident,'exit_code':0,'samples':3,'sample_max_RSS_plus_Huge_KiB':1000,'child_rusage_maxrss_units':'KiB','child_rusage_maxrss':1000}))
  (p/('samples.rank-%d.jsonl.gz'%rank)).write_bytes(b'fixture')
  text=('FACTOR_HISTORY k/teacher/step/count/coeff/ready: %d 782 12544 0 0 T\nFACTOR_SVD_QUALITY n/residual/unit_u/unit_v/smin/smax: 16 1e-15 1e-15 1e-15 0.1 1\n'%(rank+1)) if rank<80 else 'EXX_SERVICE_COMPUTE\nEXX_SERVICE_STOP\n'
  if rank==0:
   text+='end SALMON\n'
   for step in range(12537,12545):text+='FACTOR_HISTORY step/predicted: %d %s\n%d %.15e 0 0 0 32 0\n'%(step,'F' if step%16==0 else 'T',step,step*.03*.024188843265857)
  (d/('stdout.1.%d'%rank)).write_text(text)
 def ac(t):
  tw=882.447392251770225;omega=.0569614493722652368;tt=t-tw/2
  return -5.338e-9*1e6/omega*math.cos(math.pi*tt/tw)**2*math.sin(omega*tt)
 rt=[]
 for step in range(12537,12545):
  a=ac(step*.03);e=-(a-ac((step-1)*.03))/.03;rt.append([step*.03,0,0,a,0,0,e,0,0,a,0,0,e,0,0,0])
 (p/'Si_rt.data').write_text('\n'.join(' '.join(map(str,row)) for row in rt))
 (p/'Si_rt_energy.data').write_text('\n'.join('%s 0 0'%(step*.03) for step in range(12537,12545)))
 (p/'Si_nex.data').write_text('%s 0.1 0.1'%(12544*.03))
 (p/'Si_ovlp.data').write_text('12544\n'+'\n'.join(str(i)+' '+' '.join(['2']*16+['0']*16) for i in range(1,81)))
 def run():return subprocess.run([sys.executable,str(w/'post_gate.py')],cwd=p,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 assert run().returncode==0,run().stderr.decode()
 print('VALID_RESTART_GATE_ACCEPTED')
 cases=[]
 f=d/'stdout.1.159';old=f.read_text();f.unlink();assert run().returncode!=0;f.write_text(old);cases.append('missing helper')
 for fname,modify,label in [('stdout.1.1',lambda s:s+'\nNaN\n','nonfinite stdout'),('stdout.1.1',lambda s:s.replace('782 12544','782 12528'),'wrong teacher'),('stdout.1.0',lambda s:s.replace('12537 T','12537 F'),'wrong prediction')]:
  f=d/fname;old=f.read_text();f.write_text(modify(old));assert run().returncode!=0;f.write_text(old);cases.append(label)
 f=p/'Si_rt.data';old=f.read_text();s=old.splitlines();x=s[0].split();x[3]='0';s[0]=' '.join(x);f.write_text('\n'.join(s));assert run().returncode!=0;f.write_text(old);cases.append('reset field')
 f=p/'summary.rank-1.json';old=f.read_text();x=json.loads(old);x['sample_max_RSS_plus_Huge_KiB']=7*1024**2;f.write_text(json.dumps(x));assert run().returncode!=0;f.write_text(old);cases.append('rank over capacity')
 f=p/'Si_rt_energy.data';old=f.read_text();f.write_text('\n'.join(old.splitlines()[1:]));assert run().returncode!=0;f.write_text(old);cases.append('missing energy step')
 assert run().returncode==0
 print('RESTART_GATE_7_NEGATIVE_REJECTIONS_PASS',cases)
