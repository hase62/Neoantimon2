#!/usr/bin/env python3
"""Actual file:// Chrome sessions exercising explicit saved review files.

python candidate_review.py --acceptance DIR --output DIR --chrome EXECUTABLE
Requires optional Playwright in the validation environment, never Core runtime.
"""
import argparse
import hashlib
import json
from pathlib import Path
from playwright.sync_api import sync_playwright, expect


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['acceptance','output','chrome']:p.add_argument('--'+arg,required=True)
    args=p.parse_args();root=Path(args.acceptance).resolve();out=Path(args.output).resolve()
    out.mkdir(parents=True,exist_ok=False);errors=[];remote=[];checks=[];contexts=[]
    def check(name,value):
        if not value:raise AssertionError(name)
        checks.append({'check':name,'status':'PASS'})
    with sync_playwright() as pw:
        browser=pw.chromium.launch(executable_path=args.chrome,headless=True,downloads_path=str(out))
        def session(file):
            contexts.append(file)
            ctx=browser.new_context(accept_downloads=True);page=ctx.new_page()
            page.on('pageerror',lambda e:errors.append(str(e)))
            page.on('request',lambda r:remote.append(r.url) if r.url.startswith(('http:','https:')) else None)
            page.goto((root/file).as_uri());return ctx,page
        def data(page):return page.locator('#report-data').evaluate('(e)=>JSON.parse(e.textContent)')
        def select(page,id):page.locator('#rows button').filter(has_text=id).click()
        def save(page,button,name):
            with page.expect_download() as download:page.locator(button).click()
            path=out/name;download.value.save_as(path);return path
        def load(page,file,expected):
            page.locator('#review-notice').evaluate('(e)=>{e.textContent=""}')
            page.locator('#review-load').set_input_files(str(file))
            expect(page.locator('#review-notice')).to_contain_text(expected)
        ctx,page=session('private-v1.html');d=data(page);candidates=d['tables']['pmhc_candidates']
        accepted=next(c for c in candidates if c['peptide_sequence']=='ACDEFGHIK')
        held=next(c for c in candidates if c['candidate_id']!=accepted['candidate_id'])
        check('two candidates / three record Origins',len(candidates)==2 and len(d['tables']['origins'])==3)
        select(page,accepted['candidate_id'])
        check('both duplicate FASTA headers traceable',page.locator('#selected-origins').locator('.origin-record').count()==2 and
              'first original header' in page.locator('#selected-origins').inner_text() and
              'second duplicate with separate provenance' in page.locator('#selected-origins').inner_text())
        check('original score / native percentile / separate reference percentile displayed',
              'affinity_nM' in page.locator('#selected-predictions').inner_text() and
              'affinity_percentile' in page.locator('#selected-predictions').inner_text() and
              '15.625' in page.locator('#detail-body').inner_text())
        comment='=HYPERLINK("invalid") <img src=x onerror="globalThis.PWNED=1"> literal comment'
        page.locator('#review-state').select_option('accepted');page.locator('#review-comment').fill(comment)
        page.locator('#review-record').click()
        select(page,held['candidate_id']);page.locator('#review-state').select_option('held')
        page.locator('#review-comment').fill('保留: software verification');page.locator('#review-record').click()
        review=save(page,'#review-save','review-v1.json');saved=json.loads(review.read_text())
        check('JSON contains accepted / held history',len(saved['history'])==2 and {h['state']for h in saved['history']}=={'accepted','held'})
        page.locator('#review-filter').select_option('accepted');check('review-state filtering',page.locator('#rows tr').count()==1)
        table=save(page,'#review-export','accepted-v1.tsv')
        check('accepted TSV contains one confirmed candidate / excludes comments',len(table.read_text().splitlines())==2 and comment not in table.read_text())
        ctx.close()
        # Completely independent browser context: no cookies, localStorage or page state.
        ctx,page=session('private-v1.html');check('fresh session begins unreviewed','採用' not in page.locator('#rows').inner_text())
        load(page,review,'レビューを再開しました');select(page,accepted['candidate_id'])
        check('saved JSON resumes actual decision and literal comment',page.locator('#review-state').input_value()=='accepted' and page.locator('#review-comment').input_value()==comment)
        check('comment did not execute markup',page.locator('img').count()==0 and page.evaluate('globalThis.PWNED===undefined'))
        # Reject invalid files without changing already restored state.
        for name,edit in [('unknown',lambda x:x.update(schema_version='unknown')),('foreign',lambda x:x['project_identity'].update(project_run_id='different')),
                          ('corrupt',lambda x:x['history'][0].update(evidence_sha256='0'*64))]:
            value=json.loads(review.read_text());edit(value);file=out/(name+'.json');file.write_text(json.dumps(value))
            load(page,file,'読み込み拒否');check(name+' JSON rejected atomically',page.locator('#review-state').input_value()=='accepted')
        malformed=out/'malformed.json';malformed.write_text('{broken');load(page,malformed,'読み込み拒否')
        check('malformed JSON rejected atomically',page.locator('#review-comment').input_value()==comment)
        ctx.close()
        ctx,page=session('private-v2.html');load(page,review,'レビューを再開しました');select(page,accepted['candidate_id'])
        check('new analysis marks prior acceptance stale','根拠更新・再確認が必要' in page.locator('#candidate-review').inner_text())
        check('old/current evidence snapshots are both readable','判断時の根拠' in page.locator('#candidate-review').inner_text() and '現在の根拠' in page.locator('#candidate-review').inner_text())
        stale=save(page,'#review-export','accepted-stale.tsv');check('stale accepted candidate withheld from confirmed output',len(stale.read_text().splitlines())==1)
        page.locator('#review-record').click();confirmed=save(page,'#review-export','accepted-reconfirmed.tsv')
        check('explicit current-evidence confirmation enables accepted output',len(confirmed.read_text().splitlines())==2)
        final=save(page,'#review-save','review-reconfirmed.json');check('reconfirmation appends history',len(json.loads(final.read_text())['history'])==3)
        page.screenshot(path=str(out/'private-reconfirmed.png'),full_page=True);ctx.close()
        ctx,page=session('private-subset.html');load(page,review,'適用対象外');partial=save(page,'#review-save','review-partial-roundtrip.json')
        check('partial HTML preserves all out-of-scope history',json.loads(partial.read_text())['history']==saved['history'])
        select(page,data(page)['tables']['pmhc_candidates'][0]['candidate_id'])
        check('partial display does not spuriously change included-candidate evidence',
              '根拠更新・再確認が必要' not in page.locator('#candidate-review').inner_text())
        ctx.close()
        ctx,page=session('private-foreign.html');load(page,review,'読み込み拒否')
        check('separate base project cannot auto-apply review','採用' not in page.locator('#rows').inner_text());ctx.close()
        ctx,page=session('share.html');shared=data(page)
        check('sharing mode contains no review identity/history',not shared.get('review') and comment not in page.content())
        check('sharing mode disables review editing/load/export',not page.locator('#review-controls').is_visible() and page.locator('#review-filter').is_disabled())
        check('sharing mode has aliases instead of original candidate IDs',all(c['candidate_id'] not in page.content() for c in candidates))
        ctx.close()
        ctx,page=session('private-molecular.html');molecular=data(page)
        candidate=molecular['tables']['pmhc_candidates'][0];select(page,candidate['candidate_id'])
        check('molecular project still exposes authoritative linked Origins',page.locator('#selected-origins .origin-record').count()>0)
        page.locator('#review-state').select_option('held');page.locator('#review-comment').fill('Molecular Origin review retained')
        page.locator('#review-record').click();molecular_review=save(page,'#review-save','review-molecular.json');ctx.close()
        ctx,page=session('private-molecular.html');load(page,molecular_review,'レビューを再開しました');select(page,candidate['candidate_id'])
        check('same human review round trip works on existing molecular project',page.locator('#review-state').input_value()=='held')
        check('review does not mutate embedded scientific candidate records',data(page)['tables']['pmhc_candidates']==molecular['tables']['pmhc_candidates'])
        ctx.close();check('no page JavaScript errors',not errors);check('no external network requests',not remote)
        receipt={'status':'PASS','browser_version':browser.version,'scheme':'file://','isolated_contexts':len(contexts),
                 'checks':checks,'page_errors':errors,'remote_requests':remote,
                 'artifacts':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()}}
        browser.close()
    (out/'receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':len(checks)},ensure_ascii=False))


if __name__=='__main__':main()
