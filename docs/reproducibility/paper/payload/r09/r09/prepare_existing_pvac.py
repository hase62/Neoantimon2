from pathlib import Path
import argparse,subprocess,json,shutil,os,hashlib
P=Path(__file__).resolve().parent
a=argparse.ArgumentParser();a.add_argument('--python',required=True);v=a.parse_args();env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');env.pop('PYTHONPATH',None)
probe="import sys,json,pvactools,importlib.metadata as m;d=m.distribution('pvactools');print(json.dumps(dict(package=list(pvactools.__path__)[0],distribution=str(d._path),entry=sys.prefix+'/bin/pvacseq',version=m.version('pvactools'))))"
x=json.loads(subprocess.check_output([v.python,'-B','-c',probe],env=env));assert x['version']=='7.1.5';src=Path(x['package']);expected=json.loads((P/'PVACTOOLS_EXPECTED_SOURCE.json').read_text());need=sum(r['bytes']for r in expected['files']);assert shutil.disk_usage(P).free>25*1024**3+need
assert not(P/'pvac-isolated').exists()and not(P/'bin').exists(),'Use a fresh candidate directory'
for r in expected['files']:
 q=src.parent/r['file'];assert hashlib.sha256(q.read_bytes()).hexdigest()==r['sha256'],str(q)
tgt=P/'pvac-isolated';tgt.mkdir();shutil.copytree(src,tgt/'pvactools',ignore=shutil.ignore_patterns('__pycache__','*.pyc'));shutil.copytree(x['distribution'],tgt/Path(x['distribution']).name);(P/'bin').mkdir();shutil.copyfile(x['entry'],P/'bin/pvacseq');(P/'PREPARED_ENVIRONMENT.json').write_text(json.dumps(x,indent=2)+'\n');print('Prepared unchanged existing pVACtools7.1.5; no downloads/install or predictor inference')
