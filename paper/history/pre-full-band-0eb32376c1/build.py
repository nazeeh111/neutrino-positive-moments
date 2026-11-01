"""Build current research draft, retain successful output, render fresh pages."""
import hashlib,json,os,subprocess,tempfile,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
env=os.environ.copy();env['MPL_IGNORE_SYSTEM_FONTS']='1';env['MPLCONFIGDIR']=str(ROOT/'.plot-cache')
subprocess.run([__import__('sys').executable,str(ROOT/'paper/figures.py')],env=env,check=True)
import typst,pymupdf
current=ROOT/'paper/draft-current.pdf'
try:typst.compile(str(ROOT/'paper/draft.typ'),output=str(current),root=str(ROOT),ignore_system_fonts=True)
except typst.TypstError as e:print(e.diagnostic);raise
sha=hashlib.sha256(current.read_bytes()).hexdigest()
preview_root=ROOT/'paper/preview';preview_root.mkdir(exist_ok=True)
preview=Path(tempfile.mkdtemp(prefix=sha[:10]+'-',dir=preview_root))
renderer=shutil.which('pdftoppm')
if renderer is None:
    local=Path('/Users/nazeeh/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm')
    if local.exists():renderer=str(local)
if renderer:
    subprocess.run([renderer,'-png','-r','110',str(current),str(preview/'page')],check=True)
    rendering_engine='Poppler'
else:
    with pymupdf.open(current) as doc:
        for i,page in enumerate(doc):page.get_pixmap(matrix=pymupdf.Matrix(1.5,1.5)).save(preview/f'page-{i+1}.png')
    rendering_engine='PyMuPDF'
with pymupdf.open(current) as doc:
 text='\n'.join(p.get_text() for p in doc);pages=len(doc)
if pages<4 or len(text.split())<1800:raise RuntimeError('draft text/page sanity failed')
current.replace(ROOT/'paper/research-draft.pdf')
meta={'pdf_sha256':sha,'pages':pages,'word_count':len(text.split()),'preview':str(preview),'renderer':rendering_engine,
      'source_sha256':hashlib.sha256((ROOT/'paper/draft.typ').read_bytes()).hexdigest()}
(ROOT/'results/paper-build.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta))
