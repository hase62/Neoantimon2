"""Run unchanged Neoantimon2 and real pVACtools reporting CLIs on fixed controlled input."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,shutil,subprocess,sys,tarfile
P=Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now().astimezone().isoformat()
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 a=argparse.ArgumentParser();a.add_argument('--work',required=True);a.add_argument('--pvac-python',required=True);v=a.parse_args();W=Path(v.work).resolve();W.mkdir(parents=True,exist_ok=False);(W/'logs').mkdir();commands=[]
 frozen=json.loads((P/'FREEZE_BEFORE_EXECUTION.json').read_text());assert all(h(P/n)==s for n,s in frozen['files'].items())
 if(P/'SHA256SUMS.txt').exists():
  for l in(P/'SHA256SUMS.txt').read_text().splitlines():s,n=l.split('  ',1);assert h(P/n)==s
 vendor=P.parent/'pvac-isolated';entry=P.parent/'bin/pvacseq';assert vendor.is_dir()and entry.is_file()
 lib=W/'R-library';lib.mkdir();env=dict(os.environ,R_LIBS_USER=str(lib),R_MAX_VSIZE='6Gb',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONPATH=str(vendor))
 def call(label,cmd):
  if shutil.disk_usage(W).free<25*1024**3:raise RuntimeError('25 GiB free disk reserve guard')
  start=now();(W/'PROGRESS.md').write_text('RUNNING '+label+'\n'+start+'\n')
  with(W/'logs'/(label+'.stdout')).open('wb')as out,(W/'logs'/(label+'.stderr')).open('wb')as err:
   try:res=subprocess.run([str(x)for x in cmd],env=env,cwd=P,stdout=out,stderr=err,timeout=3600);code=res.returncode
   except subprocess.TimeoutExpired:code=124
  commands.append(dict(stage=label,command=[str(x)for x in cmd],started=start,finished=now(),exit_status=code,environment={'PYTHONPATH':str(vendor),'R_LIBS_USER':str(lib),'PYTHONDONTWRITEBYTECODE':'1'}));(W/'COMMANDS.json').write_text(json.dumps(commands,indent=2)+'\n');print(label+': exit '+str(code),flush=True);return code
 codes={};py=Path(v.pvac_python).absolute();r=shutil.which('R');rs=shutil.which('Rscript')
 probe='import sys,json,platform,importlib.metadata as m; import pvactools.lib.aggregate_all_epitopes as a;print(json.dumps(dict(Python=sys.version,executable=sys.executable,prefix=sys.prefix,architecture=platform.machine(),pvactools=m.version("pvactools"),imported_aggregate=a.__file__,versions={n:m.version(n)for n in ["pvactools","numpy","pandas","scipy","pysam"]}),indent=2))'
 codes['pvac_environment']=call('pvac-environment',[py,'-B','-c',probe]);assert codes['pvac_environment']==0
 ev=json.loads((W/'logs/pvac-environment.stdout').read_text());assert ev['pvactools']=='7.1.5'and Path(ev['imported_aggregate']).is_relative_to(vendor)
 for label,args in [('help-binding',['binding_filter','-h']),('help-aggregate',['generate_aggregated_report','-h'])]:codes[label]=call(label,[py,'-B',entry,*args])
 dest=W/'source';dest.mkdir()
 with tarfile.open(P/'sources/patched.tar.gz')as t:
  for m in t.getmembers():assert not Path(m.name).is_absolute()and'..'not in Path(m.name).parts and not m.issym()and not m.islnk()
  t.extractall(dest,filter='data')
 source=next(dest.iterdir());shashes={p.relative_to(source).as_posix():h(p)for p in source.rglob('*')if p.is_file()};(W/'CORE_SOURCE_HASHES.json').write_text(json.dumps(shashes,indent=2)+'\n')
 codes['R_install']=call('R-install',[r,'CMD','INSTALL','--library='+str(lib),source])
 if codes['R_install']==0:
  run=W/'neo';run.mkdir();shutil.copytree(P/'inputs',run/'inputs');(run/'runs').mkdir();(run/'observed').mkdir();codes['neo_run']=call('neo-run',[rs,'--vanilla',P/'scripts/neo_run.R',run,'R09-resolved-contexts'])
 # The pVAC INPUT table is not an observed pVAC report. Only the following real CLI outputs are counted as observations.
 pv=W/'pvac';pv.mkdir();inp=P/'inputs/R09.CONTROL_INPUT.all_epitopes.tsv'
 codes['pvac_binding']=call('pvac-binding-filter',[py,'-B',entry,'binding_filter',inp,pv/'R09.binding.filtered.tsv'])
 for mode,args in [('default',[]),('affinity-explicit',['-m2','ic50','--aggregate-inclusion-binding-threshold','5000','--aggregate-inclusion-count-limit','1'])]:
  out=pv/mode;out.mkdir();codes['pvac_aggregate_'+mode]=call('pvac-aggregate-'+mode,[py,'-B',entry,'generate_aggregated_report',inp,out/'R09.aggregated.tsv',*args])
 if codes.get('neo_run')==0:codes['neo_reopen']=call('neo-separate-reopen',[rs,'--vanilla',P/'scripts/neo_reopen.R',W/'neo'])
 for n,s in frozen['files'].items():assert h(P/n)==s,n
 assert shashes=={p.relative_to(source).as_posix():h(p)for p in source.rglob('*')if p.is_file()}
 (W/'EXECUTION_RECEIPT.json').write_text(json.dumps({'finished':now(),'exit_statuses':codes,'Core_source_unchanged':True,'frozen_inputs_expectations_unchanged':True,'actual_pvac_module':ev,'comparison_scope':'Downstream controlled-input import/reporting; no pvacseq run, VEP, trained-model inference or unknown-phase comparison','original_table_is_input':True,'scientific_query_verdict':'Pending independent normalization/query checks'},indent=2)+'\n');(W/'PROGRESS.md').write_text('COMPLETED executions; independent queries pending\n'+now()+'\n')
 sys.exit(0 if all(x==0 for x in codes.values())else 1)
if __name__=='__main__':main()
