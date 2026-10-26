#!/usr/bin/env python3
"""Offline local Markdown/anchor and editable figure checks; no analysis/network.

Run from any directory: python3 scripts/development/validate_user_docs.py
Optionally add --registry /path/to/Neoantimon2-Converters to check its public docs.
"""
from pathlib import Path
import argparse, importlib.util, json, re, struct, tempfile
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET


def slug(text):
    return re.sub(r'[^\w\- ]', '', text.lower()).replace(' ', '-')


def prose(text):
    lines=[]; fence=None
    for line in text.splitlines():
        match=re.match(r'^\s*(`{3,}|~{3,})',line)
        if match:
            if fence is None:fence=match.group(1)[0]
            elif fence==match.group(1)[0]:fence=None
            continue
        if fence is None:lines.append(line)
    return '\n'.join(lines)


def anchors(path):
    text=prose(path.read_text()); ids=set(re.findall(r'<(?:a|h[1-6])\b[^>]*\bid=["\']([^"\']+)',text));counts={}
    for value in re.findall(r'^#{1,6}\s+(.+?)\s*#*$',text,re.M):
        name=slug(re.sub(r'\[([^\]]+)\]\([^)]*\)',r'\1',value));n=counts.get(name,0);counts[name]=n+1;ids.add(name+(f'-{n}' if n else ''))
    return ids


def links(root):
    # Only version-controlled Markdown plus new documentation; exclude .git/cache.
    paths=[p for p in root.rglob('*.md') if '.git' not in p.parts and '__pycache__' not in p.parts]
    errors=[];count=0
    for path in paths:
        text=prose(path.read_text())
        targets=re.findall(r'!?\[[^\]]*\]\(([^\s()]+)(?:\s+"[^"]*")?\)',text)
        targets+=re.findall(r'<(?:a|img)\b[^>]*\b(?:href|src)=["\']([^"\']+)',text)
        for target in targets:
            url=urlsplit(target)
            if url.scheme or url.netloc:continue
            dest=(root/url.path.lstrip('/') if url.path.startswith('/') else path.parent/unquote(url.path)) if url.path else path
            count+=1
            if not dest.exists():errors.append(f'{path.relative_to(root)}: missing {target}')
            elif url.fragment and dest.suffix=='.md' and unquote(url.fragment) not in anchors(dest):errors.append(f'{path.relative_to(root)}: missing anchor {target}')
    return dict(markdown_files=len(paths),local_links=count,errors=errors)


def figure(root,base,renderer):
    source=root/'docs/assets/source'/f'{base}.drawio';svg=root/'docs/assets'/f'{base}.svg'
    doc=ET.parse(source);model=doc.find('.//mxGraphModel');cells=model.findall('./root/mxCell');ids=[c.get('id') for c in cells]
    assert len(set(ids))==len(ids)
    assert all(c.get('source') in ids and c.get('target') in ids for c in cells if c.get('edge')=='1')
    spec=importlib.util.spec_from_file_location('figure_renderer',renderer);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory(prefix='neo-doc-figure-') as temp:
        out=Path(temp)/'figure.svg';module.render(source,out);assert out.read_bytes()==svg.read_bytes(),'SVG/source mismatch'
    image=ET.parse(svg).getroot()
    assert image.find('{http://www.w3.org/2000/svg}title') is not None
    for node in image.iter():
        assert node.tag.split('}')[-1] not in ('script','image','foreignObject')
        assert not any(k.endswith('href') and (v.startswith('http') or v.startswith('//')) for k,v in node.attrib.items())
    png=(root/'docs/assets'/f'{base}.png').read_bytes();assert png[:8]==b'\x89PNG\r\n\x1a\n'
    width,height=struct.unpack('>II',png[16:24]);assert width==3*int(model.get('pageWidth')) and height==3*int(model.get('pageHeight'))
    pdf=(root/'docs/assets'/f'{base}.pdf').read_bytes();assert pdf.startswith(b'%PDF-') and b'/Font' in pdf and b'/Subtype /Image' not in pdf
    if base=='neoantimon2-input-output-examples':
        vcf=model.get('syntheticVcf');rows=vcf.splitlines();assert rows[0]=='##fileformat=VCFv4.2'
        assert len(rows[2].split('\t'))==len(rows[3].split('\t'))==10 and rows[2].split('\t')[-1]=='S1'
    return dict(name=base,editable_cells=len(cells),png=[width,height],source_svg_identical=True,vector_pdf_no_image_objects=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--registry',type=Path);args=parser.parse_args()
    core=Path(__file__).resolve().parents[2]
    result={'core':links(core),'figures':[figure(core,name,core/'scripts/development/build_workflow_figure.py') for name in ['neoantimon2-overview','neoantimon2-input-output-examples']]}
    if args.registry:
        reg=args.registry.resolve();result['registry']=links(reg);result['figures'].append(figure(reg,'converter-workflow',reg/'tools/build_workflow_figure.py'))
    print(json.dumps(result,indent=2));raise SystemExit(int(bool(result['core']['errors'] or result.get('registry',{}).get('errors'))))

if __name__=='__main__':main()
