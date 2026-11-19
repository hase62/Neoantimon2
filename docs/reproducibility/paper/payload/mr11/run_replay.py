"""Private MR11 replay entry. No network, installation of dependencies, predictor or Core modification."""
from pathlib import Path
import argparse,json,hashlib,subprocess,os,sys,shutil,tarfile,time,datetime
ROOT=Path(__file__).resolve().parent
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def now():return datetime.datetime.now().astimezone().isoformat()
def main(W):
 W=W.resolve();assert not W.exists() and not W.is_relative_to(ROOT),'Use new output OUTSIDE fixed candidate';assert shutil.disk_usage(W.parent).free>=25*1024**3
 m=json.loads((ROOT/'MANIFEST.json').read_text());files={x['path']:x for x in m['files']};actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()};assert actual==set(files)|{'MANIFEST.json','SHA256SUMS'}
 for n,x in files.items():assert H(ROOT/n)==x['sha256'] and (ROOT/n).stat().st_size==x['bytes'],n
 for l in (ROOT/'SHA256SUMS').read_text().splitlines():d,n=l.split('  ',1);assert H(ROOT/n)==d,n
 W.mkdir();(W/'logs').mkdir();(W/'observed').mkdir();library=W/'library';library.mkdir();commands=[];codes={}
 env=dict(os.environ,R_LIBS_USER=str(library),R_MAX_VSIZE='6Gb',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1');env.pop('PYTHONPATH',None)
 started=now();(W/'RUN_START.json').write_text(json.dumps({'started':started,'controller_pid':os.getpid(),'new_study':'MR11_20261116_A','candidate_manifest_sha256':H(ROOT/'MANIFEST.json'),'source_archive_sha256':H(ROOT/'source/restored20api-502-files-v001.tar.gz'),'python':sys.version,'R':shutil.which('R'),'Rscript':shutil.which('Rscript'),'environment':{k:env[k] for k in ['R_LIBS_USER','R_MAX_VSIZE','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},'scope':'new bounded11, not historical18/2004 or prior6/R09'},indent=2)+'\n')
 def save(stage,state='RUNNING'):
  (W/'COMMANDS.json').write_text(json.dumps(commands,indent=2)+'\n');(W/'STATUS.json').write_text(json.dumps({'at':now(),'state':state,'stage':stage,'controller_pid':os.getpid(),'case_exit_status':codes,'disk_free_bytes':shutil.disk_usage(W).free},indent=2)+'\n');(W/'PROGRESS.md').write_text(f'{state}: {stage}\n{now()}\nCompleted case exits: {codes}\n')
 def call(label,cmd):
  assert shutil.disk_usage(W).free>=25*1024**3,'25GiB reserve';save(label);start=now();tic=time.monotonic()
  with (W/'logs'/(label+'.stdout')).open('wb') as out,(W/'logs'/(label+'.stderr')).open('wb') as err:
   p=subprocess.Popen([str(x) for x in cmd],cwd=ROOT,env=env,stdout=out,stderr=err);r={'stage':label,'command':[str(x) for x in cmd],'cwd':str(ROOT),'pid':p.pid,'started':start};commands.append(r);save(label)
   try:code=p.wait(timeout=3600)
   except subprocess.TimeoutExpired:p.terminate();code=p.wait();r['timeout']=True
  r.update({'finished':now(),'elapsed_seconds':time.monotonic()-tic,'exit_status':code});save(label);print(label,code,flush=True);return code
 rs=shutil.which('Rscript');rb=shutil.which('R');assert rs and rb
 assert call('dependencies',[rs,'--vanilla','-e','for(p in c("jsonlite","digest","stringi","yaml","xml2"))stopifnot(requireNamespace(p,quietly=TRUE));print(sessionInfo())'])==0
 source=W/'source';source.mkdir();expected=json.loads((ROOT/'source/RESTORED_SOURCE502_HASHES.json').read_text())
 with tarfile.open(ROOT/'source/restored20api-502-files-v001.tar.gz') as t:
  for x in t.getmembers():assert not Path(x.name).is_absolute() and '..' not in Path(x.name).parts and not x.issym() and not x.islnk()
  t.extractall(source,filter='data')
 src=next(source.iterdir());assert all(H(src/n)==d for n,d in expected.items())
 assert call('install-fixed-source',[rb,'CMD','INSTALL','--library='+str(library),src])==0
 for c in json.loads((ROOT/'inputs/cases.json').read_text()):codes[c['case_id']]=call(c['case_id'],[rs,'--vanilla',ROOT/'scripts/case.R',ROOT,W,c['case_id']])
 reader=call('separate-R-reader',[rs,'--vanilla',ROOT/'scripts/reader.R',ROOT,W]);check=call('independent-saved-checker',[sys.executable,'-I',ROOT/'scripts/check_saved.py','--root',ROOT,'--out',W])
 unchanged=all(H(ROOT/n)==x['sha256'] for n,x in files.items()) and all(H(src/n)==d for n,d in expected.items());assert unchanged
 report=json.loads((W/'verification/INDEPENDENT_CHECKS.json').read_text());ok=check==0 and reader==0 and all(code==(1 if cid=='S03' else 0) for cid,code in codes.items()) and report['status']=='PASS'
 final={'status':'NEW_MOLECULAR_REPLAY_VERIFIED; HISTORICAL_STUDY18_NOT_RECOVERED' if ok else 'PARTIAL; HISTORICAL_STUDY18_NOT_RECOVERED','new_study':'MR11_20261116_A','started':started,'finished':now(),'case_exit_status':codes,'separate_reader_exit':reader,'checker_exit':check,'new_individual_checks':report['check_count'],'source_and_payload_unchanged':unchanged,'old2004_not_replayed_or_replaced':True,'previous6_R09_not_repeated':True,'new_predictor_inference':False,'source_archive':H(ROOT/'source/restored20api-502-files-v001.tar.gz'),'fixed_API_count':20,'rights':'PRIVATE_NEW_HELPERS_AND_FIXTURES_LICENSE_PENDING'}
 (W/'FINAL_EXECUTION_RECEIPT.json').write_text(json.dumps(final,indent=2)+'\n');save('finished',final['status']);print(json.dumps(final,indent=2));return 0 if ok else 1
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();raise SystemExit(main(a.out))
