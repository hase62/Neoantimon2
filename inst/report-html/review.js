'use strict';
// Explicit JSON persistence. No localStorage, cookies, network, or scientific mutation.
globalThis.NeoCandidateReview = (() => {
  const version = 'candidate-review-1.0.0';
  const states = ['unreviewed', 'accepted', 'excluded', 'held'];
  const labels = {unreviewed:'未レビュー', accepted:'採用', excluded:'除外', held:'保留'};
  const plain = v => v !== null && typeof v === 'object' && !Array.isArray(v);
  const exact = (v, keys) => plain(v) && Object.keys(v).sort().join('|') === [...keys].sort().join('|');
  const nonempty = v => typeof v === 'string' && v.length > 0;
  const hash = v => typeof v === 'string' && /^[a-f0-9]{64}$/.test(v);
  const sameIds = (a,b) => JSON.stringify([...a].sort()) === JSON.stringify([...b].sort());
  async function sha256(text) {
    const bytes = new TextEncoder().encode(text);
    const digest = await globalThis.crypto.subtle.digest('SHA-256', bytes);
    return [...new Uint8Array(digest)].map(v=>v.toString(16).padStart(2,'0')).join('');
  }
  class Store {
    constructor(payload) {
      this.enabled = !!payload;
      this.payload = payload;
      this.evidence = new Map((payload?.evidence || []).map(e=>[e.candidate_id,e]));
      this.history = [];
      this.revision = 0;
    }
    requirePrivate() { if (!this.enabled) throw new Error('レビュー操作は明示的な非公開HTMLでのみ利用できます。'); }
    latest(id) { return this.history.filter(h=>h.candidate_id===id).at(-1); }
    status(id) {
      const h=this.latest(id), e=this.evidence.get(id);
      return {state:h?.state || 'unreviewed', stale:!!h && !!e &&
        (h.evidence_sha256!==e.evidence_sha256 || !sameIds(h.displayed_external_analysis_ids,e.displayed_external_analysis_ids)), in_scope:!!e};
    }
    record(id,state,comment) {
      this.requirePrivate();
      const e=this.evidence.get(id);
      if (!e || !states.includes(state) || typeof comment!=='string' || comment.length>100000)
        throw new Error('候補、状態、コメントを確認してください。');
      const h={candidate_id:id,state,comment,recorded_at_utc:new Date().toISOString(),
        evidence_sha256:e.evidence_sha256,evidence_json:e.evidence_json,
        displayed_external_analysis_ids:[...e.displayed_external_analysis_ids]};
      this.history.push(h); this.revision++; return h;
    }
    exportText() {
      this.requirePrivate();
      return JSON.stringify({schema_version:version,project_identity:this.payload.project_identity,history:this.history},null,2);
    }
    async loadText(text) {
      this.requirePrivate();
      const revision=this.revision;
      if (typeof text!=='string' || new TextEncoder().encode(text).length>10*1024*1024) throw new Error('レビューJSONが大きすぎます。');
      const x=JSON.parse(text), identity=['project_run_id','project_analysis_id','source_manifest_sha256'];
      if (!exact(x,['schema_version','project_identity','history']) || x.schema_version!==version ||
          !exact(x.project_identity,identity) || !identity.every(k=>nonempty(x.project_identity[k])) ||
          !hash(x.project_identity.source_manifest_sha256) || !Array.isArray(x.history) || x.history.length>100000)
        throw new Error('壊れたレビューJSON、または未知のスキーマ版です。画面のレビューは変更しません。');
      if (!identity.every(k=>x.project_identity[k]===this.payload.project_identity[k]))
        throw new Error('別プロジェクトまたは異なる基礎実行のレビューJSONです。適用しません。');
      const keys=['candidate_id','state','comment','recorded_at_utc','evidence_sha256','evidence_json','displayed_external_analysis_ids'];
      for (const h of x.history) {
        if (!exact(h,keys) || !nonempty(h.candidate_id) || !states.includes(h.state) ||
            typeof h.comment!=='string' || h.comment.length>100000 || !hash(h.evidence_sha256) ||
            typeof h.evidence_json!=='string' ||
            !/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d{3})?Z$/.test(h.recorded_at_utc) ||
            !Number.isFinite(Date.parse(h.recorded_at_utc)) ||
            new Date(h.recorded_at_utc).toISOString()!==h.recorded_at_utc.replace(/(?<=:\d{2})Z$/,'.000Z') ||
            !Array.isArray(h.displayed_external_analysis_ids) || !h.displayed_external_analysis_ids.every(hash) ||
            new Set(h.displayed_external_analysis_ids).size!==h.displayed_external_analysis_ids.length)
          throw new Error('不正なレビュー履歴です。画面のレビューは変更しません。');
        const evidence=JSON.parse(h.evidence_json);
        if (!plain(evidence) || evidence.candidate?.candidate_id!==h.candidate_id ||
            await sha256(h.evidence_json)!==h.evidence_sha256)
          throw new Error('判断時の根拠ハッシュが一致しません。画面のレビューは変更しません。');
      }
      // Commit only after every record, including out-of-scope history, passes.
      if (revision!==this.revision) throw new Error('読み込み中に判断が変更されました。再度読み込んでください。');
      this.history=x.history; this.revision++;
      return this.outOfScopeCount();
    }
    outOfScopeCount() { return this.history.filter(h=>!this.evidence.has(h.candidate_id)).length; }
    acceptedTsv(candidates) {
      this.requirePrivate();
      const cell=v=>{
        let s=String(v??'').replace(/[\t\r\n]/g,' ');
        if (/^\s*[=+@-]/.test(s)) s="'"+s;
        return s;
      };
      const rows=[['candidate_id','peptide','HLA','evidence_sha256','displayed_external_analysis_ids','recorded_at_utc']];
      for (const c of candidates) {
        const s=this.status(c.candidate_id),h=this.latest(c.candidate_id);
        if (s.state==='accepted' && !s.stale) rows.push([c.candidate_id,c.peptide_sequence,c.mhc_key,
          h.evidence_sha256,h.displayed_external_analysis_ids.join(';'),h.recorded_at_utc]);
      }
      return rows.map(r=>r.map(cell).join('\t')).join('\n')+'\n';
    }
  }
  return {Store,states,labels,sha256};
})();
