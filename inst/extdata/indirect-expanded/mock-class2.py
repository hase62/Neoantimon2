#!/usr/bin/env python3
# Authored controlled SOFTWARE fixture. These values are not real predictions.
import sys,pathlib
args=sys.argv[1:]
if '-h' in args:
 print('NEOANTIMON2_TEST_MOCK\n[-BA] 0\n[-context] 0\n[-inv_all] 0\n-f -a -rdir -inptype -list');sys.exit(0)
if '-list' in args:
 print('# NetMHCIIpan version 4.3k\nDRB1_0101 DRB3_0101 DRB4_0101 DRB5_0101 HLA-DQA10101-DQB10501 HLA-DPA10202-DPB11901');sys.exit(0)
val=lambda x:args[args.index(x)+1]
root=pathlib.Path(val('-rdir'));mode=(root/'data/mode.txt').read_text().strip()
print('# NetMHCIIpan version 4.3k')
header='Pos MHC Peptide Of Core Core_Rel Inverted Identity Score_EL %Rank_EL Exp_Bind'
if mode!='el-only':header+=' Score_BA %Rank_BA Affinity(nM)'
print(header)
for i,p in enumerate(pathlib.Path(val('-f')).read_text().splitlines(),1):
 inv=int('DPA1' in val('-a') or '-inv_all' in args)
 core=(p[::-1] if inv else p)[:9]
 if mode=='badcore':core='YYYYYYYYY'
 row=f'{i} {val("-a")} {p} 0 {core} 0.750 {inv} Synthetic 0.123456 20 NA'
 if mode!='el-only':row+=' 0.700000 2.00 '+('1000.00' if mode=='boundary' else '50.00')
 print(row)
