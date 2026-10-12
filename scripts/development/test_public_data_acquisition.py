#!/usr/bin/env python3
"""Offline acquisition regressions; only newly owned temporary files are used."""
import gzip
import os
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1] / 'data/download_public_data.py'
spec = importlib.util.spec_from_file_location('acquisition', SOURCE)
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix='neo-acquisition-test-')
        self.root = Path(self.scratch.name)
        self.good = self.root / 'source.fa.gz'
        self.good.write_bytes(gzip.compress(b'>chr1\nACGTACGT\n'))

    def tearDown(self):
        self.scratch.cleanup()

    def curl(self, *args, **kwargs):
        if args[0] == 'curl':
            Path(args[args.index('--output') + 1]).write_bytes(self.good.read_bytes())
        else:
            return self.real_run(*args, **kwargs)

    def acquire(self, output):
        self.real_run = a.run
        with patch.object(a, 'run', side_effect=self.curl):
            return a.acquire('https://example.invalid/release-1/source.fa.gz', output)

    def test_truncated_prior_download_is_replaced_only_after_validation(self):
        out = self.root / 'destination.fa.gz'
        out.write_bytes(self.good.read_bytes()[:-6])
        self.acquire(out)
        self.assertEqual(out.read_bytes(), self.good.read_bytes())
        self.assertFalse(list(self.root.glob('.neo-download-*')))

    def test_failed_downloader_preserves_old_destination(self):
        out = self.root / 'destination.fa.gz'
        out.write_bytes(b'partial')
        self.real_run = a.run
        def fail(*args, **kwargs):
            if args[0] == 'curl':
                Path(args[args.index('--output') + 1]).write_bytes(b'half-new')
                raise subprocess.CalledProcessError(22, 'curl')
            return self.real_run(*args, **kwargs)
        with patch.object(a, 'run', side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                a.acquire('https://example.invalid/failed', out)
        self.assertEqual(out.read_bytes(), b'partial')

    def test_successful_http_with_invalid_body_does_not_promote(self):
        self.good.write_bytes(gzip.compress(b'<html>not fasta</html>\n'))
        out = self.root / 'destination.fa.gz'
        with self.assertRaises(ValueError): self.acquire(out)
        self.assertFalse(out.exists())

    def test_complete_download_reused_without_curl(self):
        self.real_run = a.run
        with patch.object(a, 'run', wraps=a.run) as call:
            a.acquire('https://example.invalid/never-needed', self.good)
            self.assertFalse(any(c.args[0] == 'curl' for c in call.call_args_list))

    def test_truncated_plain_fasta_reconstructed_from_valid_gzip(self):
        fasta = self.root / 'reference.fa'
        fasta.write_bytes(b'>chr1\nAC\n')
        a.decompress_fasta(self.good, fasta)
        self.assertEqual(fasta.read_bytes(), b'>chr1\nACGTACGT\n')
        self.assertEqual(Path(str(fasta) + '.fai').read_text().split('\t')[1], '8')

    def test_semantic_release_and_fasta_scope_reject_valid_wrong_inputs(self):
        annotation=self.root/'annotation.gtf.gz'
        annotation.write_bytes(gzip.compress(b'##description: GRCh38 version 49\nchr1\ttest\tgene\t1\t8\t.\t+\t.\tgene_id "synthetic";\n'))
        with self.assertRaises(ValueError):a.validate(annotation,'gencode-gtf')
        annotation.write_bytes(gzip.compress(b'##description: GRCh38 version 50 (Ensembl 116)\nchr1\ttest\tgene\t1\t8\t.\t+\t.\tgene_id "synthetic";\n'))
        a.validate(annotation,'gencode-gtf')
        a.fasta_scope(self.good,{'chr1':8},exact=True)
        with self.assertRaises(ValueError):a.fasta_scope(self.good,{'chr1':20},exact=True)

    def test_profile_selection_and_explicit_destination(self):
        with patch.object(a, 'acquire_resource') as acquire:
            a.main(['--dest', str(self.root), '--profile', 'minimum', '--resource', 'imgthla'])
            self.assertEqual([c.args[0] for c in acquire.call_args_list], ['genome','gencode','mane','imgthla'])
        with patch.object(a, 'acquire_resource') as acquire:
            a.main(['--dest', str(self.root), '--profile', 'full'])
            self.assertEqual([c.args[0] for c in acquire.call_args_list], list(a.RESOURCES))

    def panel(self, contig='22', position=20000001):
        plain = self.root / 'panel.vcf'
        plain.write_text('##fileformat=VCFv4.2\n##contig=<ID=' + contig + ',length=30000000>\n'
                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                         '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n'
                         + contig + '\t' + str(position) + '\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n')
        panel = self.root / 'panel.vcf.gz'
        a.run('bcftools','view','-Oz','-o',panel,plain)
        a.run('bcftools','index','-t',panel)
        return panel

    def test_vcf_index_pair_validated_before_promotion_and_reused(self):
        panel = self.panel()
        output = self.root/'paired/panel.vcf.gz'
        real = a.run
        def copy(*args, **kwargs):
            if args[0]=='curl':
                target=Path(args[args.index('--output')+1])
                source=Path(str(panel)+'.tbi') if target.name.endswith('.tbi') else panel
                target.write_bytes(source.read_bytes())
            else: return real(*args,**kwargs)
        with patch.object(a,'run',side_effect=copy):
            a.acquire_pair('https://example.invalid/panel.vcf.gz',output,'.tbi',a.vcf_pair)
        a.vcf_pair(output)
        with patch.object(a,'run',wraps=a.run) as calls:
            a.acquire_pair('https://example.invalid/never-needed',output,'.tbi',a.vcf_pair)
            self.assertFalse(any(c.args[0]=='curl' for c in calls.call_args_list))

    def test_vcf_pair_index_failure_does_not_promote_main_file(self):
        panel=self.panel();output=self.root/'paired/panel.vcf.gz';real=a.run
        def fail_index(*args,**kwargs):
            if args[0]=='curl':
                target=Path(args[args.index('--output')+1])
                target.write_bytes(b'bad index' if target.name.endswith('.tbi') else panel.read_bytes())
            else: return real(*args,**kwargs)
        with patch.object(a,'run',side_effect=fail_index):
            with self.assertRaises(subprocess.CalledProcessError):
                a.acquire_pair('https://example.invalid/panel.vcf.gz',output,'.tbi',a.vcf_pair)
        self.assertFalse(output.exists())
        self.assertFalse(Path(str(output)+'.tbi').exists())

    def test_1000g_contig_22_regression_and_zero_prior_subset(self):
        output = self.root / 'selected'
        output.mkdir()
        bad = output / 'NA12878.chr22.20-25Mb.GRCh38.phased.vcf.gz'
        panel = self.panel()
        a.run('bcftools','view','-r','chr22:20000000-25000000','-Oz','-o',bad,panel)
        a.run('bcftools','index','-t',bad)
        result = a.extract_1000g(panel, output)
        self.assertEqual(result, {'source_contig':'22','subset_region':'22:20000000-25000000','subset_records':1})
        self.assertIn('22\t20000001',a.run('bcftools','view','-H',bad,capture=True))

    def test_1000g_chr22_also_supported_without_renaming(self):
        result = a.extract_1000g(self.panel('chr22'), self.root/'selected')
        self.assertEqual(result['source_contig'], 'chr22')
        self.assertEqual(result['subset_records'], 1)

    def test_1000g_absent_region_fails_instead_of_false_success(self):
        with self.assertRaises(ValueError):
            a.extract_1000g(self.panel(position=100),self.root/'selected')

    def test_samtools_index_124_branch_and_truncated_bam(self):
        sam = self.root / 'reads.sam'
        sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:30\n'
                       'r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n')
        bam = self.root / 'reads.bam'
        a.run('samtools','view','-b','-o',bam,sam)
        a.index_bam(bam)
        a.bam_pair(bam)
        self.assertEqual(a.run('samtools','view','-c',bam,'chr1:1-4',capture=True).strip(),'1')
        broken = self.root / 'partial.bam'
        broken.write_bytes(bam.read_bytes()[:-28])
        with self.assertRaises(subprocess.CalledProcessError): a.validate(broken,'bam')

    def test_full_rna_shell_uses_supported_samtools_index_from_empty_runs(self):
        data=self.root/'data';data.mkdir();run_root=self.root/'runs';fake=self.root/'fake-bin';fake.mkdir()
        fasta=data/'chr1.fa';fasta.write_text('>chr1\nACGTACGTACGTACGTACGT\n')
        a.ensure_fasta(fasta)
        annotation=data/'genes.gtf.gz';annotation.write_bytes(gzip.compress(b'chr1\ttest\tgene\t1\t20\t.\t+\t.\tgene_id "synthetic";\n'))
        sra=data/'SRR37849525.sra';sra.write_bytes(b'synthetic-mock-SRA')
        vcf=data/'truth.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=20>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t4\t.\tT\tC\t.\tPASS\t.\n')
        truth=data/'truth.vcf.gz';a.run('bcftools','view','-Oz','-o',truth,vcf);a.run('bcftools','index','-t',truth)
        sam=data/'reads.sam';sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\nr1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n')
        bam=data/'provided.bam';a.run('samtools','view','-b','-o',bam,sam)
        for name,body in {
            'prefetch':'#!/bin/sh\nexit 9\n',
            'vdb-validate':'#!/bin/sh\nexit 0\n',
            'fasterq-dump':"#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\nout=Path(sys.argv[sys.argv.index('--outdir')+1])\nfor n in (1,2): (out/('SRR37849525_'+str(n)+'.fastq')).write_text('@r1\\nACGT\\n+\\nIIII\\n')\n",
            'STAR':"#!/usr/bin/env python3\nimport sys,os,shutil\nfrom pathlib import Path\na=sys.argv\nif '--runMode' in a: (Path(a[a.index('--genomeDir')+1])/'genomeParameters.txt').write_text('mock-index')\nelse: shutil.copyfile(os.environ['NEO_TEST_BAM'],a[a.index('--outFileNamePrefix')+1]+'Aligned.sortedByCoord.out.bam')\n"
        }.items():
            script=fake/name;script.write_text(body);script.chmod(0o755)
        env=dict(os.environ,PATH=str(fake)+os.pathsep+os.environ['PATH'],NEO_TEST_BAM=str(bam))
        subprocess.run(['bash',str(SOURCE.parent/'download_hcc1395_full_rna_chr1.sh'),
                        '--data-root',str(data),'--run-root',str(run_root),'--sra',str(sra),
                        '--fasta',str(fasta),'--gtf',str(annotation),'--truth',str(truth),'--threads','1'],
                       env=env,check=True,stdout=subprocess.DEVNULL)
        output=run_root/'development/HCC1395_full_rna/star_chr1_alignment/HCC1395.fullRNA.chr1.Aligned.sortedByCoord.out.bam'
        a.bam_pair(output)
        self.assertTrue(sra.exists())
        self.assertFalse((data/'public_datasets').exists())

    def test_fastq_partial_record_rejected(self):
        p = self.root / 'reads.fastq'
        p.write_text('@r1\nACGT\n+\nII\n')
        with self.assertRaises(ValueError): a.validate(p,'fastq')
        p.write_text('@r1\nACGT\n+\nIIII\n')
        a.validate(p,'fastq')

    def test_imgthla_wrong_commit_never_moves_existing_checkout(self):
        dest = self.root / 'IMGTHLA'
        dest.mkdir(); (dest/'release_version.txt').write_text(a.IMGTHLA_RELEASE)
        with patch.object(a,'run',return_value='another-commit\n') as calls:
            with self.assertRaises(ValueError): a.acquire_imgthla(dest)
            self.assertEqual(len(calls.call_args_list),1)
            self.assertNotIn('checkout',calls.call_args.args)

    def test_sra_prefetch_failure_never_promotes(self):
        with patch.object(a,'require'), patch.object(a,'run',side_effect=subprocess.CalledProcessError(1,'prefetch')):
            with self.assertRaises(subprocess.CalledProcessError): a.acquire_sra('SRR37849525',self.root/'sra')
        self.assertFalse((self.root/'sra/SRR37849525/SRR37849525.sra').exists())

    def test_sra_partial_existing_replaced_after_vdb_validation(self):
        directory = self.root/'sra'; output=directory/'SRR37849525/SRR37849525.sra'
        output.parent.mkdir(parents=True);output.write_bytes(b'partial')
        def commands(*args, **kwargs):
            if args[0]=='prefetch':
                candidate=Path(args[args.index('--output-directory')+1])/'SRR37849525/SRR37849525.sra'
                candidate.parent.mkdir();candidate.write_bytes(b'complete')
            elif Path(args[1]).read_bytes()!=b'complete':
                raise subprocess.CalledProcessError(1,'vdb-validate')
        with patch.object(a,'require'),patch.object(a,'run',side_effect=commands):
            a.acquire_sra('SRR37849525',directory)
        self.assertEqual(output.read_bytes(),b'complete')


if __name__ == '__main__':
    unittest.main(verbosity=2)
