const {test}=require('node:test');
const assert=require('node:assert/strict');
require('../../inst/report-html/review.js');
const {Store,sha256}=globalThis.NeoCandidateReview;
async function payload(ids=['c1','c2'],changed=false){
  const evidence=[];
  for(const id of ids){
    const evidence_json=JSON.stringify({candidate:{candidate_id:id},predictions:changed?[{source:'External',analysis_id:'b'.repeat(64)}]:[]});
    evidence.push({candidate_id:id,evidence_json,evidence_sha256:await sha256(evidence_json),displayed_external_analysis_ids:changed?['b'.repeat(64)]:[]});
  }
  return {schema_version:'candidate-review-1.0.0',project_identity:{project_run_id:'run1',project_analysis_id:'analysis1',source_manifest_sha256:'a'.repeat(64)},evidence};
}
test('JSON round trip, literal comments, retained out-of-scope history and re-confirmation',async()=>{
  const x=new Store(await payload());
  const comment='=1+1 <img src=x onerror=alert(1)>\tcomment';
  x.record('c1','accepted',comment);x.record('c2','held','later');
  const saved=x.exportText(),y=new Store(await payload(['c1']));
  assert.equal(await y.loadText(saved),1);assert.equal(y.latest('c1').comment,comment);
  assert.equal(JSON.parse(y.exportText()).history.length,2);
  const z=new Store(await payload(['c1'],true));await z.loadText(saved);
  assert.equal(z.status('c1').stale,true);assert.equal(z.acceptedTsv([{candidate_id:'c1'}]).trim().split('\n').length,1);
  z.record('c1','accepted','checked again');assert.equal(z.status('c1').stale,false);
  const table=z.acceptedTsv([{candidate_id:'c1',peptide_sequence:'=FORMULA()',mhc_key:'+HLA'}]);
  assert.match(table,/\t'=FORMULA\(\)\t'\+HLA/);assert.doesNotMatch(table,/checked again/);
});
test('malformed, foreign, unknown, invalid timestamps, tampered evidence fail atomically',async()=>{
  const x=new Store(await payload());x.record('c1','accepted','retain');const before=x.exportText();
  const edits=[v=>{v.schema_version='unknown'},v=>{v.project_identity.project_run_id='other'},v=>{v.history[0].state='maybe'},v=>{v.history[0].comment={html:'bad'}},v=>{v.history[0].recorded_at_utc='2026-02-31T00:00:00.000Z'},v=>{v.history[0].evidence_json+=' '},v=>{v.history.push({})}];
  for(const edit of edits){const v=JSON.parse(before);edit(v);await assert.rejects(x.loadText(JSON.stringify(v)));assert.equal(x.exportText(),before);}
  await assert.rejects(x.loadText('{broken'));assert.equal(x.exportText(),before);
  const shared=new Store(null);assert.throws(()=>shared.record('c1','accepted',''));assert.throws(()=>shared.exportText());await assert.rejects(shared.loadText(before));assert.throws(()=>shared.acceptedTsv([]));
});
