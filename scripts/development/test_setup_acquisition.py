#!/usr/bin/env python3
"""Offline tests of shared setup assets and additive acquisition policies."""
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'scripts/data/download_public_data.py'
spec = importlib.util.spec_from_file_location('setup_acquisition', SOURCE)
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
CATALOG = ROOT / 'inst/setup/resource-catalog.json'


class SetupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='neo-setup-policy-')
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(patch.stopall)
        patch.object(a, '_CATALOG', None).start()
        patch.object(a, '_CATALOG_PATH', None).start()
        patch.object(a, '_CHECK_ONLY', False).start()
        patch.object(a, '_NO_OVERWRITE', False).start()

    def test_canonical_catalog_and_packaged_utility_agree(self):
        c = a.resource_catalog()
        self.assertEqual(set(c['resources']), set(a.RESOURCES))
        self.assertEqual(c['profiles']['minimum'], list(a.PROFILES['minimum']))
        self.assertEqual(a.RESOURCE_CATALOG_SHA256, hashlib.sha256(CATALOG.read_bytes()).hexdigest())
        self.assertEqual(SOURCE.read_bytes(), (ROOT / 'inst/setup/download_public_data.py').read_bytes())
        self.assertEqual(c['resources']['imgthla']['checkout']['commit'], a.IMGTHLA_COMMIT)
        self.assertEqual(c['resources']['imgthla']['release'], a.IMGTHLA_RELEASE)

    def test_reference_role_paths_match_the_existing_config_template(self):
        # Read the existing simple template as an independent role/path inventory.
        expected = {}; role = None
        for line in (ROOT / 'config/reference.primary.example.yaml').read_text().splitlines():
            if line.startswith('  reference.'):
                role = line.strip().rstrip(':').split('.', 1)[1]
            elif role and line.startswith('    path: '):
                expected[role] = line.strip().split(': ', 1)[1]
        actual = {f['role']: f['path'] for s in a.resource_catalog()['resources'].values()
                  for f in s['files'] if 'role' in f}
        self.assertEqual(actual, expected)

    def test_no_overwrite_preserves_invalid_prior_data(self):
        source = self.root / 'source.fa.gz'
        source.write_bytes(gzip.compress(b'>chr1\nACGT\n'))
        target = self.root / 'target.fa.gz'; target.write_bytes(b'user partial data')
        with patch.object(a, '_NO_OVERWRITE', True):
            with self.assertRaisesRegex(ValueError, 'retained'):
                a.acquire(source.as_uri(), target, 'fasta')
        self.assertEqual(target.read_bytes(), b'user partial data')

    def test_no_overwrite_also_protects_partial_pairs(self):
        target = self.root / 'sample.vcf.gz'; target.write_bytes(b'prior')
        with patch.object(a, '_NO_OVERWRITE', True), patch.object(a, 'validate', side_effect=ValueError('invalid')), \
             patch.object(a, 'acquire') as download:
            with self.assertRaisesRegex(ValueError, 'retained'):
                a.acquire_pair('https://example.invalid/vcf', target, '.tbi', lambda p: None)
            download.assert_not_called()
        self.assertEqual(target.read_bytes(), b'prior')
        self.assertFalse(list(self.root.glob('.neo-pair-*')))

    def test_read_only_missing_resources_do_not_create_directories(self):
        target = self.root / 'absent' / 'sample.fa.gz'
        with patch.object(a, '_CHECK_ONLY', True):
            with self.assertRaisesRegex(ValueError, 'Missing'):
                a.acquire('https://example.invalid/file', target)
        self.assertFalse(target.parent.exists())
        dest = self.root / 'subsets'
        with patch.object(a, '_CHECK_ONLY', True), patch.object(a, 'require'):
            with self.assertRaisesRegex(ValueError, 'Missing'):
                a.extract_1000g(self.root / 'panel.vcf.gz', dest)
        self.assertFalse(dest.exists())

    def test_no_overwrite_preserves_existing_acquisition_metadata(self):
        target = self.root / 'acquisition' / 'prior.json'
        target.parent.mkdir(); target.write_text('{"original": true}\n')
        with patch.object(a, '_NO_OVERWRITE', True):
            a.atomic_json(target, {'new': True})
            a.atomic_json(target.parent / 'new.json', {'new': True})
        self.assertEqual(target.read_text(), '{"original": true}\n')
        self.assertEqual(json.loads((target.parent / 'new.json').read_text()), {'new': True})

    def test_local_mane_acquisition_and_read_only_verification(self):
        c = copy.deepcopy(a.resource_catalog())
        fixture = ROOT / 'tests/fixtures/reference'
        for i, name in enumerate(('mane.tsv', 'mane.gtf')):
            source = self.root / (name + '.gz')
            source.write_bytes(gzip.compress((fixture / name).read_bytes()))
            c['resources']['mane']['files'][i]['url'] = source.as_uri()
        destination = self.root / 'data'
        with patch.object(a, '_CATALOG', c):
            result = a.acquire_resource('mane', destination, True)
            receipt = destination / 'acquisition/mane.json'; before = receipt.read_bytes()
            facts = [(f['local_relative_path'],f['observed_sha256']) for f in result['artifacts']]
            self.assertEqual(len(facts), 2)
            for relative, sha in facts:
                self.assertEqual(hashlib.sha256((destination / relative).read_bytes()).hexdigest(), sha)
            with patch.object(a, '_CHECK_ONLY', True), patch.object(a, 'acquire_resource', wraps=a.acquire_resource):
                again = a.acquire_resource('mane', destination, True)
                self.assertEqual(result['artifacts'], again['artifacts'])
            self.assertEqual(receipt.read_bytes(), before)

    def test_standalone_help_does_not_need_a_repository_or_catalog(self):
        copied = self.root / 'download_public_data.py'; shutil.copy2(SOURCE, copied)
        result = subprocess.run([sys.executable,str(copied),'--help'],capture_output=True,text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--profile {minimum,full}', result.stdout)

    def test_standalone_companion_catalog_requires_its_exact_hash(self):
        payload = CATALOG.read_bytes()
        with patch.object(Path, 'is_file', return_value=False), \
             patch.object(a, 'require'), patch.object(a, 'run', return_value=payload.decode('utf-8')):
            self.assertEqual(set(a.resource_catalog()['resources']), set(a.RESOURCES))
        with patch.object(a, '_CATALOG', None), patch.object(Path, 'is_file', return_value=False), \
             patch.object(a, 'require'), patch.object(a, 'run', return_value=payload.decode('utf-8') + ' '):
            with self.assertRaisesRegex(ValueError, 'checksum'):
                a.resource_catalog()


if __name__ == '__main__':
    unittest.main(verbosity=2)
