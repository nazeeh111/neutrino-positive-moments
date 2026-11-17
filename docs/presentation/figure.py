"""Four retained atoms, exact coordinate/weight encoding. Not an angular scan."""
from pathlib import Path
from fractions import Fraction as F
import json,html,hashlib
P=Path(__file__).parent;case=json.loads((P/'case.json').read_text())
for mobile in [False,True]:
 w,h=(375,900)if mobile else(1080,620);s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" viewBox="0 0 {w} {h}" role="img" aria-labelledby="title desc"><title id="title">An upper bound, with its witness</title><desc id="desc">Both positive species have mass 1, zero axial flux and second number moment 1/3. The plus species places weights 1/4 at u = -1 and 3/4 at u = 1/3; the minus species reflects it. At reduced wave number 1/1000, the original longitudinal factor vanishes exactly at imaginary frequency with squared growth 1333/3000000. Horizontal position is angular cosine u; stem height is atom weight on a common zero to one scale. This is one attainable completion, not a prediction for an actual background.</desc><rect width="{w}" height="{h}" fill="#f8f5ed"/>']
 def text(x,y,t,size=16,color='#183344',anchor='start'):
  s.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" fill="{color}" text-anchor="{anchor}">{html.escape(str(t))}</text>')
 def line(x1,y1,x2,y2,color='#c9cec9',width=1):s.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>')
 pad=24 if mobile else 36;text(pad,29,'NEUTRINO POSITIVE MOMENTS / EXACT WITNESS',10)
 for j,t in enumerate(['An upper bound,','with its witness.']if mobile else['An upper bound, with its witness.']):text(pad,72+j*34,t,29 if mobile else 39)
 text(pad,161 if mobile else 115,'N = 1     Q = 1/3     K = 1/1000',15)
 x0,x1=(48,327)if mobile else(65,565);height=96
 for i,name in enumerate(['plus','minus']):
  base=(329+i*202)if mobile else(301+i*177);color='#275e67'if name=='plus'else'#aa4e32';text(x0-24,base-height-29,'μ+ / positive species'if i==0 else'μ− / positive species',14,color)
  for v in [F(0),F(1)]:
   y=base-float(v)*height;line(x0,y,x1,y,'#bfc7c4'if v==0 else'#e4e3db');text(x0-13,y+4,str(v),10,'#66716e','end')
  for u in [F(-1),F(0),F(1)]:
   x=x0+float((u+1)/2)*(x1-x0);line(x,base,x,base+5);text(x,base+22,str(u),11,'#66716e','middle')
  for atom in case['species'][name]:
   u,m=F(atom['u']),F(atom['weight']);x=x0+float((u+1)/2)*(x1-x0);y=base-float(m)*height
   s.append(f'<line data-species="{name}" data-u="{u}" data-weight="{m}" x1="{x:.10f}" x2="{x:.10f}" y1="{base}" y2="{y:.10f}" stroke="{color}" stroke-width="5"/><circle cx="{x}" cy="{y}" r="4" fill="{color}"/>');text(x,y-13,str(m),14,color,'middle')
   if u not in [-1,0,1]:text(x,base+22,str(u),11,color,'middle')
  text((x0+x1)/2,base+45,'angular cosine u',11,'#66716e','middle')
 y=634 if mobile else 206;x=24 if mobile else 664
 text(x,y,'SAME SUPPLIED MOMENTS',10,'#66716e');text(x,y+31,'(mass, flux, second) = (1, 0, 1/3)',14)
 text(x,y+82,'EXACT ATTAINING MODE',10,'#66716e');text(x,y+116,'γ² = 1333 / 3000000',23)
 text(x,y+145,'D_L = 0, checked in the original factor.',12)
 text(x,y+179,'Positive weights. Exact rational arithmetic.',11)
 fy=850 if mobile else 565;text(pad,fy,'One compatible completion. Not universal instability.',11 if mobile else 14);text(pad,fy+23,'Reduced, dimensionless units; no physical wavelength assigned.',10 if mobile else 13)
 s.append('</svg>');(P/('witness-mobile.svg'if mobile else'witness.svg')).write_text(''.join(s))
(P/'figure-manifest.json').write_text(json.dumps({'source':'case.json','sourceSHA256':hashlib.sha256((P/'case.json').read_bytes()).hexdigest(),'encoding':'u position = x0+(u+1)/2*(x1-x0); stem height =96*weight; common u[-1,1] and weight[0,1]. Four actual atoms, no interpolation or scan. Equal circle radii are decoration, not mass encoding.','variants':{'witness.svg':{'x0':65,'x1':565,'bases':[301,478]},'witness-mobile.svg':{'x0':48,'x1':327,'bases':[329,531]}},'claim':'Attainability check only. Written universal upper-bound proof remains separate.','units':'Dimensionless reduced frequency and wave number; number-weighted angular moments.'},indent=2)+'\n')
