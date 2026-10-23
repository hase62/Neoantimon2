#!/usr/bin/env python3
"""Render the project's editable, uncompressed diagrams.net vector source.

Uses only Python's standard library. Supported diagram primitives: groups,
rectangles, text cells, and straight/polyline edges with normalized ports.
PNG/PDF exports are separate optional headless-Chrome operations.
"""
from pathlib import Path
import argparse
import html
import json
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

SVG='http://www.w3.org/2000/svg'

def style_map(value):
    return dict(p.split('=',1) for p in value.split(';') if '=' in p)

def render(source,destination):
    tree=ET.parse(source);model=tree.find('.//mxGraphModel')
    if model is None:raise ValueError('Expected an uncompressed diagrams.net model')
    cells={c.get('id'):c for c in model.findall('./root/mxCell')}
    if len(cells)!=len(model.findall('./root/mxCell')):raise ValueError('Duplicate cell IDs')
    width=float(model.get('pageWidth'));height=float(model.get('pageHeight'))
    def geometry(cid):
        c=cells[cid];g=c.find('mxGeometry')
        if g is None:return (0.,0.,0.,0.)
        x=float(g.get('x','0'));y=float(g.get('y','0'))
        parent=c.get('parent')
        if parent not in (None,'0','1'):
            px,py,_,_=geometry(parent);x+=px;y+=py
        return x,y,float(g.get('width','0')),float(g.get('height','0'))
    ET.register_namespace('',SVG)
    svg=ET.Element('{%s}svg'%SVG,dict(width=f'{width:g}',height=f'{height:g}',viewBox=f'0 0 {width:g} {height:g}',role='img',**{'aria-labelledby':'figure-title figure-description'}))
    ET.SubElement(svg,'title',id='figure-title').text='Neoantimon2: variants to traceable peptide–HLA candidates'
    ET.SubElement(svg,'desc',id='figure-description').text=('User-declared samples, patients and HLA plus explicit references support variant events, multiple sequence and translation paths, shared peptide sequences and peptide–HLA candidates. Separately acquired references and predictors are not bundled. Optional configured local prediction and supplied observations support native prioritization and explanations. A completed project can export queries and all origins for a separately run external model, optional conversion, strict matching and an independent derived analysis. Native results are preserved. Longitudinal, vaccine-set and restricted indirect-HLA/allo views have separate opt-in requirements. Missing evidence is neither zero nor an evaluated negative. Research outputs do not establish clinical validation.')
    defs=ET.SubElement(svg,'defs')
    for color,name in [('#176487','native'),('#a34a09','external'),('#52616d','optional')]:
        m=ET.SubElement(defs,'marker',id='arrow-'+name,markerWidth='9',markerHeight='9',refX='8',refY='4.5',orient='auto',markerUnits='userSpaceOnUse')
        ET.SubElement(m,'path',d='M 0 0 L 8 4.5 L 0 9 Z',fill=color)
    ET.SubElement(svg,'rect',width=f'{width:g}',height=f'{height:g}',fill='#ffffff')
    # Containers precede connectors; opaque labels/cards are painted last.
    ordered=[c for c in cells.values() if c.get('vertex')=='1' and style_map(c.get('style','')).get('neoLayer')=='background']
    ordered += [c for c in cells.values() if c.get('edge')=='1']
    ordered += [c for c in cells.values() if c.get('vertex')=='1' and style_map(c.get('style','')).get('neoLayer')!='background']
    for c in ordered:
        cid=c.get('id');s=style_map(c.get('style',''))
        if s.get('shape')=='group' or 'group' in c.get('style','').split(';'):continue
        group=ET.SubElement(svg,'g',id=cid)
        if c.get('edge')=='1':
            def port(node,prefix):
                x,y,w,h=geometry(node);return x+float(s.get(prefix+'X','0.5'))*w,y+float(s.get(prefix+'Y','0.5'))*h
            start=port(c.get('source'),'exit');end=port(c.get('target'),'entry')
            points=[start]+[(float(p.get('x')),float(p.get('y'))) for p in c.findall('mxGeometry/Array[@as="points"]/mxPoint')]+[end]
            color=s.get('strokeColor','#52616d');marker='native' if color=='#176487' else 'external' if color=='#a34a09' else 'optional'
            attrs=dict(points=' '.join(f'{x:g},{y:g}' for x,y in points),fill='none',stroke=color,**{'stroke-width':s.get('strokeWidth','2.3'),'stroke-linejoin':'round','stroke-linecap':'round','marker-end':'url(#arrow-'+marker+')'})
            if s.get('dashed')=='1':attrs['stroke-dasharray']='7 5'
            ET.SubElement(group,'polyline',attrs);continue
        x,y,w,h=geometry(cid);shape=s.get('shape','rectangle')
        if shape not in ('text','rectangle'):raise ValueError('Unsupported shape: '+shape)
        if shape=='rectangle':
            attrs=dict(x=f'{x:g}',y=f'{y:g}',width=f'{w:g}',height=f'{h:g}',fill=s.get('fillColor','#ffffff'),stroke=s.get('strokeColor','#c5cfd7'),**{'stroke-width':s.get('strokeWidth','1.4')})
            if s.get('rounded')=='1':attrs.update(rx='10',ry='10')
            if s.get('dashed')=='1':attrs['stroke-dasharray']='7 5'
            ET.SubElement(group,'rect',attrs)
        value=c.get('value','')
        if not value:continue
        size=float(s.get('fontSize','18'));align=s.get('align','left')
        anchor={'left':'start','center':'middle','right':'end'}[align]
        tx=x if align=='left' else x+w/2 if align=='center' else x+w
        attrs=dict(x=f'{tx:g}',y=f'{y+size*.94:g}',fill=s.get('fontColor','#233746'),**{'font-family':'Arial, Helvetica, sans-serif','font-size':f'{size:g}','text-anchor':anchor,'data-box':f'{x:g} {y:g} {w:g} {h:g}'})
        font_style=int(s.get('fontStyle','0'))
        if font_style&1:attrs['font-weight']='700'
        if font_style&2:attrs['font-style']='italic'
        text=ET.SubElement(group,'text',attrs)
        for i,line in enumerate(value.splitlines()):
            ET.SubElement(text,'tspan',x=f'{tx:g}',y=f'{y+size*.94+i*size*1.2:g}').text=line
    ET.indent(svg,space='  ');destination.parent.mkdir(parents=True,exist_ok=True)
    ET.ElementTree(svg).write(destination,encoding='utf-8',xml_declaration=True)

def browser_exports(svg_path, png, pdf, browser, work_dir):
    """Render exports and retain inspectable previews; Chrome is optional."""
    if work_dir is None:
        raise ValueError('--work-dir is required for browser exports/previews')
    work_dir = work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    executable = browser or shutil.which('chromium') or shutil.which('google-chrome')
    if not executable:
        candidate = Path('/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
        executable = str(candidate) if candidate.is_file() else None
    if not executable:
        raise ValueError('Pass --browser with a Chrome/Chromium executable')
    width, height = (float(x) for x in ET.parse(svg_path).getroot().get('viewBox').split()[2:])
    svg = svg_path.read_text(encoding='utf-8').split('?>', 1)[-1]
    metrics_script = """<script>window.onload=()=>{
      let issues=[];for(let t of document.querySelectorAll('text[data-box]')){
        let [x,y,w,h]=t.dataset.box.split(' ').map(Number),b=t.getBBox();
        if(b.x<x-1||b.y<y-3||b.x+b.width>x+w+1||b.y+b.height>y+h+1)
          issues.push({id:t.parentElement.id,box:[x,y,w,h],actual:[b.x,b.y,b.width,b.height]});
      }let p=document.createElement('pre');p.id='metrics';
      p.textContent=JSON.stringify({issues,textElements:document.querySelectorAll('text').length});
      document.body.append(p);};</script>"""
    page = work_dir / 'figure.html'
    page.write_text('<!doctype html><meta charset="utf-8"><title>Neoantimon2 workflow overview</title><style>@page{size:'+
      f'{width:g}px {height:g}px'+
      ';margin:0}html,body{margin:0;background:white}svg{display:block;width:100%;height:auto}pre{display:none}</style>'+
      svg+metrics_script, encoding='utf-8')
    args = [str(executable), '--headless', '--disable-gpu', '--no-first-run',
      '--no-default-browser-check', '--hide-scrollbars',
      '--user-data-dir='+str(work_dir/'chrome-profile')]
    commands = []
    def call(extra):
        cmd = args + extra + [page.as_uri()]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
        commands.append(dict(command=cmd, exit_code=result.returncode))
        (work_dir/'render-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
        if result.returncode:
            raise RuntimeError('Browser export failed: '+result.stderr[-2000:])
        return result.stdout
    dom = call(['--dump-dom'])
    (work_dir/'figure-dom.html').write_text(dom, encoding='utf-8')
    match = re.search(r'<pre id="metrics">(.*?)</pre>',dom)
    if not match:
        raise RuntimeError('Browser did not produce SVG text measurements')
    metrics = json.loads(html.unescape(match.group(1)))
    (work_dir/'text-metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    if metrics['issues']:
        raise ValueError('Text exceeds source boxes; inspect text-metrics.json')
    for pixels in (1200,900,600):
        target = work_dir/f'overview-{pixels}.png'
        call([f'--window-size={pixels},{round(pixels*height/width)}', '--screenshot='+str(target)])
        if not target.is_file():
            raise RuntimeError('Browser produced no preview image')
    if png:
        png = png.resolve(); png.parent.mkdir(parents=True,exist_ok=True)
        call([f'--window-size={round(width*3)},{round(height*3)}','--screenshot='+str(png)])
        if not png.is_file():
            raise RuntimeError('Browser produced no PNG export')
    if pdf:
        pdf = pdf.resolve(); pdf.parent.mkdir(parents=True,exist_ok=True)
        call(['--no-pdf-header-footer','--print-to-pdf='+str(pdf)])
        if not pdf.is_file():
            raise RuntimeError('Browser produced no PDF export')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--png',type=Path,help='Optional 3x-resolution PNG export')
    p.add_argument('--pdf',type=Path,help='Optional vector PDF export')
    p.add_argument('--browser',help='Chrome/Chromium executable (for exports only)')
    p.add_argument('--work-dir',type=Path,help='Retain browser previews and measurements here')
    args=p.parse_args()
    render(args.source,args.output)
    if args.png or args.pdf or args.work_dir:
        browser_exports(args.output,args.png,args.pdf,args.browser,args.work_dir)
