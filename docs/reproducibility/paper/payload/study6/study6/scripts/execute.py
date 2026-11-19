"""Bounded closure execution: isolated source/library, no downloads or prediction engines."""
from pathlib import Path
import sys,json,subprocess,os,shutil,hashlib,datetime,tarfile,argparse
P=Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now().astimezone().isoformat()
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 a=argparse.ArgumentParser();a.add_argument('--work',required=True);a.add_argument('--case-ids',nargs='*');v=a.parse_args();W=Path(v.work).resolve();W.mkdir(parents=True,exist_ok=False);(W/'logs').mkdir();cmds=[]
 freeze=json.loads((P/('FREEZE_BEFORE_FRESH_V3.json'if(P/'FREEZE_BEFORE_FRESH_V3.json').exists()else'FREEZE_BEFORE_CORE_V2.json'if(P/'FREEZE_BEFORE_CORE_V2.json').exists()else'FREEZE_BEFORE_CORE.json')).read_text())
 for line in (P/'SHA256SUMS.txt').read_text().splitlines():
  sha,n=line.split('  ',1);assert h(P/n)==sha,n
 for n,sha in freeze['files'].items():assert h(P/n)==sha,n
 lib=W/'library';lib.mkdir();env=dict(os.environ,R_LIBS_USER=str(lib),R_MAX_VSIZE='6Gb',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 (W/'EXECUTION_ENV.json').write_text(json.dumps({'started':now(),'source_archive_sha256':h(P/'sources/patched.tar.gz'),'source_identity':json.loads((P/'SOURCE_IDENTITY.json').read_text()),'R':shutil.which('R'),'Rscript':shutil.which('Rscript'),'python':sys.version,'platform':sys.platform,'env_overrides':{n:env[n]for n in ['R_LIBS_USER','R_MAX_VSIZE','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},'limits':'25 GiB free reserve; one sequential R child;6GiB vector heap is not total RSS cap'},indent=2)+'\n')
 def call(stage,cmd):
  if shutil.disk_usage(W).free<25*1024**3:raise RuntimeError('25 GiB free reserve guard')
  (W/'PROGRESS.md').write_text('RUNNING '+stage+'\n'+now()+'\n')
  started=now()
  with(W/'logs'/(stage+'.stdout')).open('wb')as out,(W/'logs'/(stage+'.stderr')).open('wb')as err:
   try:r=subprocess.run(cmd,env=env,cwd=P,stdout=out,stderr=err,timeout=3600);code=r.returncode
   except subprocess.TimeoutExpired:code=124
  rec={'stage':stage,'command':cmd,'started':started,'finished':now(),'exit_status':code,'env_R_LIBS_USER':str(lib)};cmds.append(rec);(W/'COMMANDS.json').write_text(json.dumps(cmds,indent=2)+'\n');print(stage+': exit '+str(code),flush=True);return code
 rs=shutil.which('Rscript');r=shutil.which('R')
 assert call('dependencies',[rs,'--vanilla','-e','for(p in c("jsonlite","digest","stringi","yaml","xml2"))stopifnot(requireNamespace(p,quietly=TRUE));print(sessionInfo())'])==0
 dest=W/'source';dest.mkdir()
 with tarfile.open(P/'sources/patched.tar.gz')as t:
  for m in t.getmembers():assert not Path(m.name).is_absolute()and'..'not in Path(m.name).parts and not m.issym()and not m.islnk()
  t.extractall(dest,filter='data')
 source=next(dest.iterdir());sourcehash={f.relative_to(source).as_posix():h(f)for f in source.rglob('*')if f.is_file()};(W/'EXTRACTED_SOURCE_HASHES.json').write_text(json.dumps(sourcehash,indent=2)+'\n')
 assert call('install',[r,'CMD','INSTALL','--library='+str(lib),str(source)])==0
 run=W/'run';run.mkdir();shutil.copytree(P/'inputs',run/'inputs');(run/'runs').mkdir();(run/'observed').mkdir()
 codes={}
 for c in [x for x in json.loads((P/'inputs/cases.json').read_text())if not v.case_ids or x['case_id']in v.case_ids]:codes[c['case_id']]=call(c['case_id'],[rs,'--vanilla',str(P/'scripts/run_case.R'),str(run),c['case_id']])
 reopen=call('separate-process-reopen',[rs,'--vanilla',str(P/'scripts/reopen_saved.R'),str(run)])
 for n,sha in freeze['files'].items():assert h(P/n)==sha,n
 for n,sha in sourcehash.items():assert h(source/n)==sha,n
 summary={'finished':now(),'case_exit_status':codes,'reopen_exit_status':reopen,'source_and_frozen_inputs_unchanged':True,'qualification':'Execution exit status is separate from scientific expected-value comparison; see independent review outputs'};(W/'EXECUTION_RECEIPT.json').write_text(json.dumps(summary,indent=2)+'\n');(W/'PROGRESS.md').write_text('COMPLETED execution; inspect scientific comparison separately\n'+now()+'\n')
 review_code=call('independent-saved-review',[sys.executable,str(P/'scripts/review_saved.py'),str(W)])if(P/'scripts/review_saved.py').exists()else None
 summary['independent_saved_review_exit_status']=review_code;summary['execution_and_review_ok']=all(x==0 for x in codes.values())and reopen==0 and review_code==0
 (W/'EXECUTION_RECEIPT.json').write_text(json.dumps(summary,indent=2)+'\n')
 sys.exit(0 if summary['execution_and_review_ok']else 1)
if __name__=='__main__':main()
