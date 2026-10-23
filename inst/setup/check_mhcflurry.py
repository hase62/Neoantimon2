"""Read-only runtime probe; no downloads, installs, environment or model changes."""
import importlib
import json
import sys

result = dict(python=sys.executable, python_version=sys.version.split()[0])
try:
    try:
        from importlib import metadata
        version = metadata.version('mhcflurry')
        requires = metadata.requires('mhcflurry') or []
    except ImportError:
        # Python 3.7 installations can still run historical MHCflurry releases.
        import pkg_resources
        distribution = pkg_resources.get_distribution('mhcflurry')
        version = distribution.version
        requires = [str(r) for r in distribution.requires()]
    result['version'] = version
    # Older MHCflurry omits TensorFlow from install_requires because multiple
    # distributions provide it. Merely passing pip's metadata checks is insufficient.
    backend = 'torch' if any(r.lower().startswith('torch') for r in requires) or (
        tuple(int(x) for x in version.split('.')[:2]) >= (2, 3)) else 'tensorflow'
    result['backend'] = backend
    for module in ('mhcflurry', 'numpy', 'pandas', backend):
        importlib.import_module(module)
    result['success'] = True
except Exception as error:
    result.update(success=False, error=type(error).__name__ + ': ' + str(error))
print('NEO_SETUP_JSON:' + json.dumps(result))
sys.exit(0 if result['success'] else 1)
