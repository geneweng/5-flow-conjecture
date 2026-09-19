import re, subprocess, math, html, sys, os

def build_page(md_path, out_path, eyebrow, nav_html, standalone=True):
  md = open(md_path).read()
  # drop the H1 and italic byline and the plain Contents list; we render those ourselves
  title = re.match(r'^# (.*)', md).group(1)
  md = re.sub(r'^# .*?\n', '', md, count=1)
  bm = re.search(r'^\*((?:Prepared|Working notes).*?)\*\s*$', md, re.M)
  byline = bm.group(1) if bm else ''
  md = re.sub(r'^\*(?:Prepared|Working notes).*?\*\s*$', '', md, count=1, flags=re.M)
  md = re.sub(r'## Contents\n\n(?:\d+\..*\n)+', '', md)
  body = subprocess.run(['pandoc','--from','markdown','--to','html5','--mathjax','--wrap=none'],
                        input=md, capture_output=True, text=True, check=True).stdout
  # table of contents from h2s
  toc = re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body)
  toc_html = '<nav class="toc" aria-label="Contents"><p class="eyebrow">Contents</p><ol>' + ''.join(
      f'<li><a href="#{i}">{re.sub(r"^\d+\.\s*","",t)}</a></li>' for i,t in toc) + '</ol></nav>'
  # wrap tables for horizontal scroll
  body = body.replace('<table>', '<div class="tbl"><table>').replace('</table>', '</table></div>')
  # the conjecture block: first blockquote-like bold statement
  body = body.replace('<p><strong>Conjecture (Tutte, 1954).</strong> Every bridgeless graph admits a nowhere-zero 5-flow.</p>',
    '<p class="conj"><span class="lbl">Conjecture (Tutte, 1954)</span> Every bridgeless graph admits a nowhere-zero 5-flow.</p>')
  # Petersen graph SVG
  R1, R2 = 40, 18
  outer = [(60+R1*math.sin(2*math.pi*i/5), 62-R1*math.cos(2*math.pi*i/5)) for i in range(5)]
  inner = [(60+R2*math.sin(2*math.pi*i/5), 62-R2*math.cos(2*math.pi*i/5)) for i in range(5)]
  edges = [(outer[i],outer[(i+1)%5]) for i in range(5)] + [(outer[i],inner[i]) for i in range(5)] + [(inner[i],inner[(i+2)%5]) for i in range(5)]
  svg = '<svg class="petersen" viewBox="0 0 120 120" width="120" height="120" role="img" aria-label="The Petersen graph"><g stroke="var(--accent)" stroke-width="1.6" fill="none">'
  svg += ''.join(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}"/>' for a,b in edges)
  svg += '</g><g fill="var(--ink)">' + ''.join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3"/>' for x,y in outer+inner) + '</g></svg>'
  head = '<title>' + html.escape(title) + '</title>' + '''
  <link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&family=Source+Sans+3:wght@400;600;700&display=swap">
  <script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)'],['$','$']],displayMath:[['\\\\[','\\\\]'],['$$','$$']]},svg:{fontCache:'global'}};</script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.js"></script>
  <style>
  :root{--paper:#f6f6f3;--paper2:#ecede8;--ink:#1c1f23;--ink2:#4f5760;--rule:#cfd3cc;--accent:#2b5c85;--accent2:#9a4b2f;--code:#e9ebe6;}
  @media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--paper:#15181c;--paper2:#1e2227;--ink:#e7e5df;--ink2:#a2a9b1;--rule:#353b42;--accent:#86b6dc;--accent2:#d99570;--code:#22272d;}}
  :root[data-theme="dark"]{--paper:#15181c;--paper2:#1e2227;--ink:#e7e5df;--ink2:#a2a9b1;--rule:#353b42;--accent:#86b6dc;--accent2:#d99570;--code:#22272d;}
  html{color-scheme:light dark}
  body{background:var(--paper);color:var(--ink);font-family:"Source Serif 4",Georgia,"Times New Roman",serif;font-size:17px;line-height:1.55;padding-block:0 4rem;padding-inline:clamp(16px,4vw,40px);margin:0}
  .wrap{max-width:46rem;margin:0 auto}
  header.masthead{display:grid;grid-template-columns:1fr auto;gap:1.5rem;align-items:end;padding-block:3rem 1.5rem;border-bottom:2px solid var(--ink)}
  header.masthead h1{font-family:"Source Sans 3","Helvetica Neue",Arial,sans-serif;font-weight:700;font-size:clamp(1.9rem,4.5vw,2.8rem);line-height:1.05;letter-spacing:-.01em;margin:0 0 .6rem;text-wrap:balance}
  header.masthead .byline{color:var(--ink2);font-size:.95rem;margin:0;max-width:34rem}
  .eyebrow{font-family:"Source Sans 3",Arial,sans-serif;font-weight:600;font-size:.74rem;letter-spacing:.12em;text-transform:uppercase;color:var(--ink2);margin:0 0 .4rem}
  .petersen{flex:none}
  .sitenav{font-family:"Source Sans 3",Arial,sans-serif;font-size:.9rem;color:var(--ink2);margin:1rem 0 0}
.conj{font-size:1.2rem;line-height:1.4;margin:2rem 0;padding:1rem 1.2rem;border-left:3px solid var(--accent);background:var(--paper2)}
  .conj .lbl{font-family:"Source Sans 3",Arial,sans-serif;font-weight:600;font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--accent);display:block;margin-bottom:.3rem}
  nav.toc{margin:2rem 0 2.5rem;padding:1rem 1.2rem;border:1px solid var(--rule)}
  nav.toc ol{margin:0;padding-left:1.4rem;columns:2;column-gap:2rem;font-family:"Source Sans 3",Arial,sans-serif;font-size:.95rem;line-height:1.5}
  nav.toc li{break-inside:avoid;padding:.1rem 0}
  @media (max-width:520px){nav.toc ol{columns:1}header.masthead{grid-template-columns:1fr}.petersen{width:84px;height:84px}}
  h2{font-family:"Source Sans 3",Arial,sans-serif;font-weight:700;font-size:1.45rem;line-height:1.2;margin:3rem 0 1rem;padding-top:1rem;border-top:1px solid var(--rule);text-wrap:balance}
  h3{font-family:"Source Sans 3",Arial,sans-serif;font-weight:600;font-size:1.08rem;margin:2rem 0 .6rem;text-wrap:balance}
  p{margin:0 0 1rem}
  a{color:var(--accent);text-decoration-thickness:1px;text-underline-offset:2px}
  a:focus-visible{outline:2px solid var(--accent2);outline-offset:2px}
  blockquote{margin:1.2rem 0;padding:.7rem 1.1rem;border-left:3px solid var(--accent);background:var(--paper2)}
  blockquote p:last-child{margin-bottom:0}
  ul,ol{padding-left:1.4rem}li{margin-bottom:.35rem}
  hr{border:0;height:0;margin:0}
  .tbl{overflow-x:auto;margin:1.2rem 0}
  table{border-collapse:collapse;width:100%;font-family:"Source Sans 3",Arial,sans-serif;font-size:.92rem;line-height:1.4;font-variant-numeric:tabular-nums}
  th{text-align:left;font-weight:600;font-size:.78rem;letter-spacing:.08em;text-transform:uppercase;color:var(--ink2);border-bottom:2px solid var(--ink);padding:.45rem .6rem .35rem}
  td{border-bottom:1px solid var(--rule);padding:.45rem .6rem;vertical-align:top}
  td:first-child{font-weight:600}
  code{font-family:"SF Mono",Menlo,Consolas,monospace;font-size:.86em;background:var(--code);padding:.05em .3em;border-radius:2px}
  strong{font-weight:600}
  mjx-container[display="true"]{margin:.6rem 0;overflow-x:auto;overflow-y:hidden}
  section.refs ul{padding-left:1.2rem;font-size:.93rem}
  @media (prefers-reduced-motion:reduce){*{scroll-behavior:auto}}
  </style>
  '''
  page = head + f'''<div class="wrap">
  <header class="masthead"><div><p class="eyebrow">{html.escape(eyebrow)}</p><h1>{html.escape(title)}</h1><p class="byline">{html.escape(byline)}</p></div>{svg}</header>
  {toc_html}
  {body}
  </div>'''
  if standalone:
    page = '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n' + page.replace('<div class="wrap">', '</head>\n<body>\n<div class="wrap">', 1) + '\n</body>\n</html>\n'
  open(out_path, 'w').write(page)
  print(out_path, len(page), 'bytes;', len(toc), 'sections')

if __name__ == '__main__':
  root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  nav = '<nav class="sitenav"><a href="index.html">Survey</a> · <a href="mod5-orientation.html">Route: modulo 5-orientations</a> · <a href="oddness.html">Route: oddness 6</a> · <a href="c6-model-spec.html">(C6) model</a> · Experiments <a href="exp-B-defect.html">B</a> <a href="exp-C-rounding.html">C</a> <a href="exp-D-counting.html">D</a> · <a href="https://github.com/geneweng/5-flow-conjecture">GitHub</a></nav>'
  build_page(f'{root}/survey/5-flow-survey.md', f'{root}/docs/index.html', 'Survey · September 2026', nav)
  build_page(f'{root}/notes/mod5-orientation.md', f'{root}/docs/mod5-orientation.html', 'Working notes · September 2026', nav)
  build_page(f'{root}/notes/oddness.md', f'{root}/docs/oddness.html', 'Working notes · September 2026', nav)
  for extra in ('c6-model-spec', 'exp-B-defect', 'exp-C-rounding', 'exp-D-counting', 'c6-independent-check'):
    if os.path.exists(f'{root}/notes/{extra}.md'): build_page(f'{root}/notes/{extra}.md', f'{root}/docs/{extra}.html', 'Working notes · September 2026', nav)
  build_page(f'{root}/survey/5-flow-survey.md', f'{root}/survey/5-flow-survey.html', 'Survey · September 2026', '', standalone=False)
