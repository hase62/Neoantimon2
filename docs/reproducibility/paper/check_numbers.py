#!/usr/bin/env python3
"""Reaggregate frozen observations and saved receipts. Prints JSON; no study execution."""
from collections import Counter, defaultdict
from contextlib import contextmanager
from decimal import Decimal
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import re
import zipfile

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'observations'
SNAPSHOT = 'f4fb8eb9dff5af88b6cb73fc19aa582ef415753f45df3851cc60d555c527ef23'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_json(name):
    return json.loads((DATA / name).read_text())


@contextmanager
def tsv(name):
    path = DATA / name
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        yield csv.DictReader(handle, delimiter='\t')


def records(name):
    with tsv(name) as stream:
        return list(stream)


def unique(rows, key):
    values = [row[key] for row in rows]
    require(all(values) and len(values) == len(set(values)), f'Non-unique or empty {key}')
    return set(values)


def verify():
    expected = json.loads((ROOT / 'observations-sha256.json').read_text())
    observed = {p.relative_to(DATA).as_posix() for p in DATA.rglob('*') if p.is_file()}
    require(observed == expected.keys(), 'Saved observation file inventory differs')
    for name, item in expected.items():
        path = DATA / name
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        require(digest == item['sha256'] and path.stat().st_size == item['bytes'],
                f'Saved observation identity mismatch: {name}')


def oracle_case(row):
    expected, observed = row['expected'], row['observed']
    if row['group'] == 'fusion':
        wanted = set(map(tuple, expected['windows']))
        literal = (wanted == set(map(tuple, observed['windows']))
                   and bool(observed['contexts'])
                   and all(v['protein'] == expected['protein'] and v['dna'] == expected['dna']
                           for v in observed['contexts']))
        return literal and expect_executed(row), not wanted
    proteins = {x['mutant_protein'] for x in expected if x['mutant_protein'] is not None}
    got_proteins = {x['protein'] for x in observed if x['protein'] is not None}
    wanted = {(l, p) for x in expected for l, _, p in x['windows']}
    got = {(l, p) for y in observed for l, _, p, _, _ in y['windows']}
    literal = proteins == got_proteins and wanted == got
    for x in expected:
        if x['mutant_protein'] is None:
            continue
        context = next((y for y in observed if y['protein'] == x['mutant_protein']
                        and ('transcript_id' not in x or y['transcript'] == x['transcript_id'])), None)
        if context is None:
            continue  # Historical literal comparison used the first compatible context if available.
        actual = {(l, start, p) for l, start, p, _, _ in context['windows']}
        literal &= context['dna'] == x['mutant_transcript'] and actual == set(map(tuple, x['windows']))
        normal = x['normal_protein']
        if len(normal) == len(x['mutant_protein']):
            literal &= all(wt is None or wt == normal[start-1:end]
                           for _, start, _, end, wt in context['windows'])
    return literal and expect_executed(row), not wanted


def expect_executed(row):
    return row['status'] == 'EXECUTED'


def oracle():
    cases = read_json('oracle/cases.json')
    unique(cases, 'case_id')
    result = Counter()
    qualified = []
    for row in cases:
        match, empty = oracle_case(row)
        result[row['group']] += 1
        result['literal_matches'] += match
        result['zero_windows_expected'] += empty
        result['zero_windows_matched'] += empty and match
        if not match:
            qualified.append(row['case_id'])
    graphs = {}
    sets = {}
    for case in ('same_codon_cis','same_codon_trans'):
        graph = read_json(f'oracle/{case}-pipeline.json')
        origins = unique(graph['origins'],'origin_id')
        unique(graph['candidates'],'candidate_id')
        sets[case] = {x['peptide_sequence'] for x in graph['origins']}
        links = [(r['candidate_id'],o) for r in graph['candidates'] for o in r['origin_ids']]
        require(len(links)==len(set(links)) and all(o in origins for _,o in links),'Oracle graph link mismatch')
        graphs[case] = dict(sequences=len(sets[case]),Origins=len(origins),links=len(links))
    graphs['sequence_intersection'] = len(set.intersection(*sets.values()))
    case_by_id = {row['case_id']: row for row in cases}
    def observed_sequences(case):
        return {window[2] for context in case_by_id[case]['observed'] for window in context['windows']}

    baseline = observed_sequences('single_missense')
    germline = observed_sequences('germline_somatic_cis')
    public_germline = read_json('oracle/public-background.json')
    background = {case: dict(Origins=len(rows), sequences=len({r['peptide_sequence'] for r in rows}),
                             same_as_single_missense={r['peptide_sequence'] for r in rows} == baseline,
                             phase_states=Counter(r['phase_status'] for r in rows))
                  for case, rows in public_germline.items()}
    return dict(cases=len(cases), **result, qualified=qualified, archived_cis_trans=graphs,
                low_level_cis_background=dict(shared=len(baseline & germline),
                                              gained=len(germline - baseline), lost=len(baseline - germline)),
                public_germline_background=background, synthetic_isoform=isoform())


def isoform():
    saved = read_json('oracle/long-read.json')
    genome = fasta(DATA / 'oracle/long-read-genome.fasta')
    codons = dict(zip(map(''.join, itertools.product('TCAG', repeat=3)),
                      'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
    proteins = {}
    for form in ('wt', 'mt'):
        record = saved[form]
        structure = record['structure']
        contig = genome[structure['contig']]
        dna = ''.join(contig[e['start']-1:e['end']] for e in structure['exons'])
        require(dna == record['dna'], f'Isoform exon projection differs for {form}')
        amino_acids = (codons[dna[i:i+3]] for i in range(0, len(dna)-2, 3))
        proteins[form] = ''.join(itertools.takewhile(lambda aa: aa != '*', amino_acids))
    require(set(saved['mt']['proteins']) == {proteins['mt']}, 'Isoform translation differs')
    candidate, = saved['candidates']
    origins = unique(saved['origins'], 'origin_id')
    links = {(row['candidate_id'], row['origin_id']) for row in saved['candidate_origin_links']}
    require(set(json.loads(candidate['origin_ids'])) == origins
            and links == {(candidate['candidate_id'], origin) for origin in origins}
            and len(links) == len(saved['candidate_origin_links'])
            and candidate['peptide_sequence'] == proteins['mt'], 'Incomplete isoform linkage')
    require({o for group in saved['source_Origins'].values() for o in group} == origins
            and sum(map(len, saved['source_Origins'].values())) == len(origins)
            and {row['peptide_sequence'] for row in saved['origins']} == {proteins['mt']},
            'Isoform source Origins differ')
    return dict(proteins=proteins, source_domains={k: len(v) for k, v in saved['source_Origins'].items()},
                Origins=len(origins), pMHC=len(saved['candidates']), links=len(links))


def fasta(path):
    saved = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:]
            require(name not in saved, 'Duplicate FASTA header')
            saved[name] = ''
        else:
            require(name is not None, 'FASTA without a header')
            saved[name] += line.strip()
    return saved


def vep_mutant(annotation):
    normal = annotation['WildtypeProtein']
    if annotation['FrameshiftSequence']:
        return annotation['FrameshiftSequence'], 'literal_vep_frameshift'
    if 'frameshift_variant' in annotation['Consequence']:
        return None, 'unavailable'
    position = annotation['Protein_position'].split('/')[0].split('-')
    start, end = int(position[0])-1, int(position[-1])
    ref, alt = annotation['Amino_acids'].split('/')
    ref = '' if ref == '-' else ref
    alt = '' if alt == '-' else alt
    require(normal[start:end] == ref, 'VEP amino acid interval does not match saved WT')
    method = 'derived_stop_gain' if 'stop_gained' in annotation['Consequence'] else 'derived_edit'
    return (normal[:start]+alt+normal[end:]).split('*',1)[0].split('X',1)[0], method


def pvac():
    native = read_json('pvac/native.json')
    vep = read_json('pvac/vep-matched.json')
    native_keys = {(r['event_id'], r['transcript_id']): r for r in native}
    locus = {(r['contig'], r['pos'] if len(r['ref'])>len(r['alt']) else r['pos']-1,
              r['pos']-1+len(r['ref']), r['ref'], r['alt']):r['event_id'] for r in native}
    annotations = {(r['event_id'],r['transcript']):r['annotations'] for r in vep}
    require(len(native_keys)==len(native) and len(annotations)==len(vep), 'Duplicate context')
    require(set(native_keys)==set(annotations) and all(len(v)==1 for v in annotations.values()),
            'VEP match must be exact and unique')
    converted = {}
    for r in records('pvac/converted.tsv'):
        event = locus[r['chromosome_name'],int(r['start']),int(r['stop']),r['reference'],r['variant']]
        key = event,r['transcript_name']
        require(key not in converted and key in native_keys, 'Invalid converted context')
        converted[key] = r['index']
    fragments = defaultdict(list)
    for header, sequence in fasta(DATA/'pvac/protein-fragments.fasta').items():
        if header.startswith('MT.'):
            fragments[header.removeprefix('MT.')].append(sequence)
    result = Counter()
    for k, index in converted.items():
        a = annotations[k][0]
        n = native_keys[k]
        require(a['Feature']==n['transcript_id'], 'Transcript version mismatch')
        mt, mode = vep_mutant(a)
        require(mt is not None, 'Converted context lacks a full comparison view')
        result[mode] += 1
        result['full_mutant_equal'] += mt == n['tumor_protein']
        result['full_wt_equal'] += a['WildtypeProtein'] == n['patient_normal_protein']
        result['literal_mt_fragment_contexts'] += bool(fragments[index])
        result['fragment_substring_contexts'] += bool(fragments[index]) and all(f in mt for f in fragments[index])
    excluded = set(native_keys)-set(converted)
    require(all('splice' in annotations[k][0]['Consequence'] and vep_mutant(annotations[k][0])[0] is None
                for k in excluded), 'An excluded context lacks the recorded unavailability reason')
    events = records('pvac/selected-events.tsv')
    unique(events,'event_id')
    return dict(source_variants=len(events), variant_classes=Counter(r['category'] for r in events),
                Neo_contexts=len(native), VEP_matched=len(vep),
                upstream_vep_rows=len(records('pvac/vep-unfiltered.tsv')),
                shared=len(converted), shared_variants=len({k[0] for k in converted}),
                unavailable=len(excluded), unavailable_variants=len({k[0] for k in excluded}), **result)


def endpoints(text, limit):
    values = [int(v) for v in text.split(',')] if text else []
    require(len(set(values))==len(values) and all(0<=v<limit for v in values),
            'Duplicate or unwired endpoint index')
    return len(values)


def full_graph():
    events=records('full-graph/events.tsv.gz')
    unique(events,'event_id')
    require(len({(r['input_id'],r['line'],r['selected_alt_index']) for r in events})==len(events),
            'Source record keys are not unique')
    origins=records('full-graph/origins.tsv.gz')
    unique(origins,'origin_id')
    origin_events=sum(endpoints(r['event_row_indices'],len(events)) for r in origins)
    candidates=records('full-graph/candidates.tsv.gz')
    unique(candidates,'candidate_id')
    links=sum(endpoints(r['origin_row_indices'],len(origins)) for r in candidates)
    return dict(Events=len(events), Origins=len(origins), pMHC=len(candidates), links=links,
                Origin_Event_links=origin_events, source_records_by_input=Counter(r['input_id'] for r in events),
                original_driver=read_json('full-graph/driver-exit.json'),
                original_controller=read_json('full-graph/controller-exit.json'))


def graph(folder, suffix=''):
    def get(name):
        return records(f'{folder}/{name}.tsv{suffix}')
    events, origins, candidates = get('events'),get('origins'),get('candidates')
    oids=unique(origins,'origin_id')
    cids=unique(candidates,'candidate_id')
    eids=unique(events,'event_id')
    if origins and 'source_event_ids' in origins[0]:
        require(all(set(json.loads(row['source_event_ids'])) <= eids for row in origins),
                'Pilot Origin references an unknown Event')
    link_rows=get('candidate_origin_links')
    edges={(r['candidate_id'],r['origin_id']) for r in link_rows}
    require(len(edges)==len(link_rows)
            and all(c in cids and o in oids for c,o in edges),'Invalid candidate–Origin links')
    sequences={r['peptide_sequence'] for r in origins}
    counts=dict(Events=len(events),Origins=len(origins),sequences=len(sequences),pMHC=len(candidates),links=len(edges))
    return counts,sequences,candidates


def s18_origin_status(row):
    dna = json.loads(row['origin_dna_evidence'])
    phase = json.loads(row['origin_phase_condition'])['status']
    require(isinstance(dna, list) and dna, 'S18 Origin has no required-Event DNA evidence')
    unique(dna, 'event_id')
    require(all(type(event[k]) is int and event[k] >= 0
                for event in dna for k in ('tumor_alt', 'normal_alt')),
            'S18 ALT count is missing or invalid')
    require(phase in {'RESOLVED', 'CONDITIONAL'}, 'S18 phase condition is missing or invalid')
    if any(event['tumor_alt'] < 5 or event['normal_alt'] > 1 for event in dna):
        return 'ineligible'
    return 'not_evaluated' if phase == 'CONDITIONAL' else 'eligible'


def s18():
    old = {r['case_id']:r for r in read_json('saved/study18-native.json') if r['version']=='patched'}
    study6 = {(r['case_id'], r['candidate_id'], r['representative_origin_id']):
              (r['single_Origin_external_query_status'],r['all_Origins_native_status'])
              for r in records('saved/study6-D_all_representatives.tsv')}
    native6 = {(r['case_id'],r['candidate_id']):
               (r['HLA'],r['actual_native_status'],r['actual_selected'].lower())
               for r in records('saved/study6-native_target_decisions.tsv')}
    rows = records('saved/Neoantimon2_S18_observations.tsv')
    cases = set(old) | {case for case, _, _ in study6}
    require({r['case_id'] for r in rows} == cases and len(cases) == 6, 'S18 condition coverage differs')
    require({r['sequence'] for r in rows} == {'MKACDFGH'}, 'S18 target peptide differs')
    groups = defaultdict(list)
    for row in rows:
        expected_scope = 'study18_patched_saved' if row['case_id'] in old else 'study6_original_saved'
        require(row['evidence_scope'] == expected_scope, 'S18 evidence scope differs')
        groups[row['case_id'], row['HLA']].append(row)
    hlas = {'HLA-A*02:01', 'HLA-B*07:02'}
    require(set(groups) == set(itertools.product(cases, hlas)), 'S18 case/HLA grid differs')
    for case in sorted(cases):
        first, second = (groups[case,hla] for hla in sorted(hlas))
        require({r['candidate_id'] for r in first}.isdisjoint(r['candidate_id'] for r in second)
                and {r['origin_id'] for r in first} == {r['origin_id'] for r in second},
                f'S18 candidate identities or shared Origins differ: {case}')
    compared6 = set()
    status_counts = Counter()
    changed = 0
    selected = Counter()
    results = []
    for (case, hla), group in sorted(groups.items()):
        require(len(group) == 3, f'S18 Origin coverage differs: {case}/{hla}')
        unique(group, 'origin_id')
        candidate, = {r['candidate_id'] for r in group}
        native_status, = {r['observed_candidate_status'] for r in group}
        native_selected, = {r['observed_candidate_selected'] for r in group}
        require(candidate and native_selected in {'true', 'false'}, 'Invalid S18 candidate record')
        observed_origin_statuses = {r['origin_id']:r['observed_origin_status'] for r in group}
        if case in old:
            archived = {r['origin_id']:r['status'] for r in json.loads(old[case]['origins'])}
            require(archived == observed_origin_statuses and native_status == old[case]['eligibility'],
                    f'S18 legacy support differs: {case}/{hla}')
        else:
            require(native6[case,candidate] == (hla,native_status,native_selected),
                    f'S18 additional-study candidate differs: {case}/{hla}')
        derived = {r['origin_id']:s18_origin_status(r) for r in group}
        require(derived == observed_origin_statuses, f'S18 Origin status differs: {case}/{hla}')
        statuses = set(derived.values())
        aggregate = ('eligible' if 'eligible' in statuses else
                     'ineligible' if statuses == {'ineligible'} else 'not_evaluated')
        require(aggregate == native_status, f'S18 full-information status differs: {case}/{hla}')
        selected[case] += native_selected == 'true'
        for origin, status in derived.items():
            if case not in old:
                key = case,candidate,origin
                require(study6[key] == (status,native_status), f'S18 archived query differs: {key}')
                compared6.add(key)
            status_counts[status] += 1
            changed += status != native_status
        results.append(dict(case_id=case,HLA=hla,candidate_id=candidate,full_status=native_status,
                            native_selected=native_selected=='true',single_Origin_status=derived))
    require(compared6 == set(study6) and all(selected[c]==int(r['selected']) for c,r in old.items()),
            'S18 archived comparison coverage or native selection count differs')
    return dict(conditions=len(cases),HLA_candidates_per_condition=len(hlas),
                candidate_condition_groups=len(groups),representative_queries=sum(status_counts.values()),
                queries_from_legacy_conditions=sum(len(group) for (case,_),group in groups.items() if case in old),
                archived_study6_queries_matched=len(compared6),single_Origin_status_counts=status_counts,
                same_as_full=sum(status_counts.values())-changed,different_from_full=changed,groups=results)


def saved():
    result={}
    for name in ('study18','integrated18'):
        cases=read_json(f'saved/{name}-checks.json')
        require(len({(r['version'],r['case_id']) for r in cases})==len(cases),'Duplicate case receipt')
        result[name]=dict(cases=len(cases), types=len({r['case_id'] for r in cases}),
                          exit_histogram=Counter(r['exit'] for r in cases),
                          checks=sum(len(r['checks']) for r in cases),
                          passed=sum(p for r in cases for _,p in r['checks']))
    old = read_json('saved/study18-checks.json')
    labels = [k for row in old for k,ok in row['checks'] if ok]
    raw_count = sum(k.startswith('raw_evidence_identity_') for k in labels)
    omitted = sum(k.startswith('omitted_is_missing') for k in labels)
    result['study18_controlled_import'] = dict(raw_rows=raw_count,omitted=omitted,
                                             candidate_entries=raw_count+omitted)
    result['integrated18']['baseline_reused'] = len(read_json('saved/integration-baseline-reuse.json'))
    result['study18_first_four']= {r['case_id']:json.loads(r['counts'])
        for r in read_json('saved/study18-executions.json')
        if r['version']=='patched' and r['case_id'].startswith('NV-EE-')}
    result['study18_native']=[dict(case=r['case_id'],eligibility=r['eligibility'],selected=int(r['selected']),
                                    per_HLA_Origin_status=Counter(o['status'] for o in json.loads(r['origins'])))
                              for r in read_json('saved/study18-native.json') if r['version']=='patched']
    result['study18_readers']={name:Counter(r['status'] for r in read_json(f'saved/old-{name}-readers.json'))
                               for name in ('baseline','patched')}
    for name in ('mr11','study6','r09-first','r09-fresh'):
        checks=read_json(f'saved/{name}-checks.json')
        result[name]=dict(checks=len(checks),passed=sum(x['PASS'] for x in checks))
    mr11_cases=read_json('saved/mr11-case-records.json')
    frozen = json.loads((ROOT / 'payload/mr11/expected/independent_expected.json').read_text())
    exits = read_json('saved/mr11-exits.json')
    require(set(frozen) == set(exits) == {row['case_id'] for row in mr11_cases}, 'MR11 case inventory differs')
    require(all(exits[case] == row['expected_exit'] for case, row in frozen.items()), 'MR11 exit contract differs')
    result['mr11']['molecular_case_classes'] = Counter(
        'invalid_annotation_rejection' if row['expected_exit'] else
        'unresolved_contract' if row.get('expected_unresolved') else 'sequence_expectation'
        for row in frozen.values())
    result['mr11']['saved_checker_classes'] = Counter(row['classification'] for row in mr11_cases)
    result['mr11']['case_Origin_and_pMHC_observed']={x['case_id']:[x.get('actual_Origins'),x.get('actual_pMHC')] for x in mr11_cases}
    result['mr11'].update(exit_histogram=Counter(read_json('saved/mr11-exits.json').values()),
                         reader_cases=[x['case_id'] for x in read_json('saved/mr11-readers.json')
                                       if x['R_process_pid']!=x['first_R_process_pid']])
    d=read_json('saved/study6-exits.json')
    result['study6'].update(exit_histogram=Counter(d['case_exit_status'].values()),reader_exit=d['reopen_exit_status'])
    native=records('saved/study6-native_target_decisions.tsv')
    origins=defaultdict(Counter)
    for r in records('saved/study6-Origin_support_decisions.tsv'):
        origins[r['case_id'],r['candidate_id']][r['observed_Origin_status']]+=1
    for case, group in itertools.groupby(sorted(native,key=lambda x:x['case_id']),key=lambda x:x['case_id']):
        all_rows=list(group)
        print_rows=Counter(x['actual_native_status'] for x in all_rows)
        result['study6'].setdefault('target_cases',{})[case]={
          'retained':len(all_rows),'selected':sum(x['actual_selected']=='True' for x in all_rows),
          'status':print_rows,'per_candidate_Origin_status':[origins[case,x['candidate_id']] for x in all_rows]}
    reps=records('saved/study6-D_all_representatives.tsv')
    result['study6']['representative_queries']=len(reps)
    result['study6']['representative_groups']=len({(x['case_id'],x['candidate_id']) for x in reps})
    result['study6']['matched_graphs']={
        c:graph(f'study6-graph/{c}')[0] for c in ('CLOSURE-A-mixed','CLOSURE-A-conditional-only')}
    result['study6']['boundary_control_proteins']=[{k:x[k] for k in ('transcript','observed_protein','first_stop_nt')}
        for x in read_json('saved/study6-molecular_contexts.json') if x['case_id']=='CLOSURE-C-valid-boundary-del']
    result['study6']['boundary_Origin_windows']=Counter(x['transcript']
        for x in read_json('saved/study6-all_Origin_windows.json') if x['case_id']=='CLOSURE-C-valid-boundary-del')
    def query_value(value):
        return json.loads(value) if value.startswith(('{','[')) else value
    first_queries = {}
    for stage in ('first','fresh'):
        q=records(f'saved/r09-{stage}-queries.tsv')
        result['r09-'+stage].update(query_rows=len(q),query_types=len({r['query_id'] for r in q}),
            routes=len({r['tool'] for r in q}))
        comparisons = {(r['tool'],r['query_id']):query_value(r['observed'])
                       for r in q if r['query_id']!='Q5'}
        matched = sum(query_value(r['observed'])==query_value(r['expected'])
                      for r in q if r['query_id']!='Q5')
        if stage=='first':
            first_queries = comparisons
        else:
            require(comparisons==first_queries,'Fresh process returned different query values')
            matched += sum(r['observed']=='same' and r['expected']=='same as first process'
                           for r in q if r['query_id']=='Q5')
        result['r09-'+stage]['matched']=matched
        relations=records(f'saved/r09-{stage}-relations.tsv')
        by_route = defaultdict(list)
        for row in relations:
            by_route[row['dataset']].append(row)
        natural_fields = [column for column in relations[0] if column != 'dataset']
        relation_sets = [{tuple(r[key] for key in natural_fields) for r in rows} for rows in by_route.values()]
        require(all(len(keys) == len(rows) for keys, rows in zip(relation_sets, by_route.values()))
                and all(keys == relation_sets[0] for keys in relation_sets), 'R09 route relations differ')
        for route, rows in by_route.items():
            for query, dropped in [('Q3a', {'chr1'}), ('Q3b', {'chrB'}), ('Q4', {'chr1', 'chrB'})]:
                subset = [r for r in rows if r['contig'] not in dropped]
                expected = dict(Events=sorted({':'.join(r[k] for k in ('contig','pos','ref','alt')) for r in subset}),
                                contexts=sorted({r['canonical_context_id'] for r in subset}),
                                pmhc=sorted({r['peptide_sequence']+'|'+r['HLA'] for r in subset}),
                                relations=len(subset), unique_sequences=sorted({r['peptide_sequence'] for r in subset}))
                require(comparisons[route, query] == expected, f'R09 recomputed query differs: {stage}/{route}/{query}')
        result['r09-'+stage]['relations_by_route']=Counter(r['dataset'] for r in relations)
        result['r09-'+stage].update(contexts=len({r['canonical_context_id'] for r in relations}),
            peptides=len({r['peptide_sequence'] for r in relations}),
            Events=len({tuple(r[k] for k in ('contig','pos','ref','alt')) for r in relations}),
            HLA_scores=sorted({(r['HLA'],r['score']) for r in relations}))
    result['r09_recorded_command_exits'] = Counter(read_json('saved/r09-exits.json').values())
    restored = read_json('saved/restored-protected-files.json')
    require(restored['before'] == restored['after'], 'Restored study preservation ledger differs')
    result['restored_protected_file_matches'] = len(restored['before'])
    result['r09_report_rows'] = {n:len(records(f'saved/r09-{n}.tsv'))
        for n in ('binding','aggregate-default','aggregate-affinity-explicit')}
    result['distribution_manifest_entries']=len(read_json('saved/distribution-sha256.json'))
    result['protected_distribution_rows']=sum(x['original_sha256']==x['distribution_sha256']
        for x in read_json('saved/distribution-input-bindings.json'))
    protected=read_json('saved/development-protected-files.json')
    require(protected['before']==protected['after'],'Development protected source changed')
    source=read_json('saved/development-source-identity.json')
    result['development_source']={**source,'protected_file_matches':len(protected['before'])}
    result['tests']={name:read_json(f'saved/test-{name}.json')['counts']
        for name in ('focused','installed-all','source-external-regression')}
    for name in ('test-legacy','test-development','r-check-integration','r-check-development'):
        result['tests'][name]=(DATA/f'saved/{name}.txt').read_text().strip()
    return result


def authentic_predictions(external):
    observed = {}
    for label,provider,seq,hla,metrics in [
        ('ba','mhcflurry-ba','peptide','allele',{'affinity_nM':'mhcflurry_affinity','affinity_percentile':'mhcflurry_affinity_percentile'}),
        ('el','bigmhc-el','pep','mhc',{'score':'BigMHC_EL'}),
        ('im','bigmhc-im','pep','mhc',{'score':'BigMHC_IM'}),
    ]:
        path=DATA/f'authentic/{label}-raw.csv'
        fingerprint=hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open(newline='') as f:
            for i,row in enumerate(csv.DictReader(f),2):
                for metric,col in metrics.items():
                    key=(provider,metric,row[seq],row[hla])
                    require(key not in observed,'Duplicate authentic raw measurement')
                    observed[key]=(row[col],f'sha256:{fingerprint}:line:{i}:column:{col}')
    fields=('provider','metric','peptide_sequence','mhc_key')
    standard = set()
    for row in itertools.chain(*(records(f'authentic/{x}-standard.tsv') for x in ('ba','el','im'))):
        key=tuple(row[k] for k in fields)
        require(key not in standard and observed[key]==(row['score'],row['source_row']),
                'Raw/standard numeric or source mismatch')
        standard.add(key)
    require(standard==observed.keys(),'Raw and standard scopes differ')
    imported={tuple(row[k] for k in fields):row for row in external}
    require(len(imported)==len(external) and imported.keys()==observed.keys(),'Raw/observed prediction coverage differs')
    require(all(imported[key]['original_score_text']==value[0] for key,value in observed.items()),
            'Saved imported original score string changed')
    return dict(measurements=len(observed), raw_standard_and_import_match=True)


def snapshot(archive, selected, pilot_queries):
    with archive.open('rb') as stream:
        require(hashlib.file_digest(stream,'sha256').hexdigest()==SNAPSHOT,'Public asset SHA-256 mismatch')
    with zipfile.ZipFile(archive) as z:
        def get(n):
            handle=io.StringIO(z.read(f'inputs/report/{n}.tsv').decode('utf-8'))
            return list(csv.DictReader(handle,delimiter='\t'))
        candidates,origins,preds=get('pmhc_candidates'),get('origins'),get('pmhc_predictions')
        edges=get('candidate_origin_links')
        causal=get('origin_event_links')
        events=get('events')
    candidate_ids = unique(candidates, 'candidate_id')
    origin_ids = unique(origins, 'origin_id')
    event_ids = unique(events, 'event_id')
    saved_source = records('pilot/source-vcf-rows.tsv')
    allele_columns = ('contig', 'pos', 'ref', 'alt')
    event_keys = {tuple(row[key] for key in allele_columns) for row in events}
    require(len(saved_source) == len(events) == len(event_keys)
            and event_keys == {tuple(row[key] for key in allele_columns) for row in saved_source},
            'Pilot source VCF allele identity differs')
    pmhc_queries = {(row['peptide_sequence'], row['mhc_key']) for row in candidates}
    require(len(pmhc_queries) == len(candidates)
            and pmhc_queries == {pair for pair in pilot_queries if pair[0] in selected},
            'Snapshot peptide-HLA membership differs from full pilot')
    require(len(edges) == len({(r['candidate_id'], r['origin_id']) for r in edges})
            and all(r['candidate_id'] in candidate_ids and r['origin_id'] in origin_ids for r in edges)
            and all(r['origin_id'] in origin_ids and r['event_id'] in event_ids for r in causal),
            'Snapshot linkage differs')
    scopes=Counter((r['record_kind'],r['provider'],r['metric']) for r in preds)
    external=[r for r in preds if r['record_kind']=='measurement' and r['source']=='External']
    trace=authentic_predictions(external)
    require(len(causal)==len({(r['origin_id'],r['event_id'],r['link_role']) for r in causal}),
            'Repeated source relation')
    kinds={'BA':('mhcflurry-ba','affinity_nM'),'EL':('bigmhc-el','score'),'IM':('bigmhc-im','score')}
    rank={}
    evaluated={}
    for label,(provider,metric) in kinds.items():
        found=[r for r in external if r['provider']==provider and r['metric']==metric and r['prediction_status']=='evaluated']
        require(bool(found),f'Missing {label} evaluated records')
        evaluated[label]={r['candidate_id'] for r in found}
        directions={r['score_direction'] for r in found}
        require(directions in ({'lower_is_better'}, {'higher_is_better'}) and
                len({(r['predictor_version'], r['run_id'], r['task']) for r in found}) == 1,
                f'Invalid {label} score scope')
        require(len(found) == len(evaluated[label]) and evaluated[label] == candidate_ids,
                f'Invalid {label} candidate coverage')
        direction=1 if directions=={'lower_is_better'} else -1
        proxy={}
        for r in sorted(found,key=lambda x:(direction*Decimal(x['raw_score']),x['mhc_key'],x['candidate_id'])):
            proxy.setdefault(r['peptide_sequence'],r)
        rank[label]=[r['peptide_sequence'] for r in sorted(proxy.values(),
             key=lambda x:(direction*Decimal(x['raw_score']),x['peptide_sequence'],x['mhc_key']))]
    overlap={}
    for n in (5,10):
        pairs={}
        for a,b in itertools.combinations(rank,2):
            x,y=set(rank[a][:n]),set(rank[b][:n])
            pairs[f'{a}/{b}']={'intersection':len(x&y),'Jaccard':len(x&y)/len(x|y)}
        overlap[n]=pairs
    seq_origins=defaultdict(set)
    for r in origins:
        seq_origins[r['peptide_sequence']].add(r['origin_id'])
    agrn='AAPCRVLPQG'
    oids=seq_origins[agrn]
    cids={r['candidate_id'] for r in candidates if r['peptide_sequence']==agrn}
    def events_for(ids):
        return {r['event_id'] for r in causal if r['origin_id'] in ids and r['link_role']=='generation_source_set_member'}
    return dict(sequences=len(seq_origins),Origins=len(origins),pMHC=len(candidates),links=len(edges),
                external=len(external),source_vcf_alleles_linked=len(event_keys),authentic_trace=trace,common_evaluated=len(set.intersection(*evaluated.values())),measurements=[[*k,v] for k,v in scopes.items()],event_records=len(events),
                linked_events=len(events_for({r['origin_id'] for r in origins})),
                origin_histogram=Counter(len(v) for v in seq_origins.values()),
                lexical_membership_matches=(selected==set(seq_origins)),overlap=overlap,
                sequences_by_gene={g:len({r['peptide_sequence'] for r in origins if r['gene_symbol']==g})
                                   for g in sorted({r['gene_symbol'] for r in origins})},
                phase=Counter(r['phase_status'] for r in origins),
                missing_Origin_gene_tpm=sum(r['gene_tpm']=='' for r in origins),
                missing_Origin_transcript_tpm=sum(r['transcript_tpm']=='' for r in origins),
                self_lookup_status=Counter(r['self_match_status'] for r in candidates),
                self_exact_match=Counter(r['self_exact_match'] for r in candidates),
                technical_eligibility=Counter(r['technical_eligibility'] for r in candidates),
                default_include=Counter(r['default_include'] for r in candidates),
                AGRN=dict(origins=len(oids),pMHC=len(cids),
                    links=sum(r['candidate_id'] in cids and r['origin_id'] in oids for r in edges),
                    measurements=sum(r['candidate_id'] in cids for r in external),events=len(events_for(oids)),
                    transcript_versions=sorted(r['transcript_id']+'.'+r['transcript_version']
                                               for r in origins if r['origin_id'] in oids)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path,help='Existing reproducibility-v1 asset')
    args=parser.parse_args()
    verify()
    result={'oracle':oracle(),'pvacseq':pvac(),'full_graph':full_graph(),'saved_replays':saved(),'s18':s18()}
    counts,seq,candidates=graph('pilot','.gz')
    pool={r['peptide_sequence'] for r in candidates if 8<=len(r['peptide_sequence'])<=12}
    require(pool=={s for s in seq if 8<=len(s)<=12},'Pilot candidate sequence pool differs')
    selected={s for n in range(8,13) for s in sorted(x for x in pool if len(x)==n)[:3]}
    result['pilot']={**counts,'Class_I_length_sequences':len(pool),'length_15_sequences':sum(len(s)==15 for s in seq),
                     'lexical_selected_sequences':sorted(selected),'Class_I_presenters':len({r['mhc_key'] for r in candidates})}
    if args.snapshot:
        pilot_queries = {(r['peptide_sequence'], r['mhc_key']) for r in candidates}
        result['snapshot']=snapshot(args.snapshot,selected,pilot_queries)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
