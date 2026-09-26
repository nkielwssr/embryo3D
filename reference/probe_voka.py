"""probe_voka.py — diagnostic du viewer : iframe interne (boutons, bulle Astuces, canvas), chargement, rotation pas à pas.
Usage : python probe_voka.py [step_px] [n_steps] [--headed] [--slow ms]"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capture_voka as C

args = [a for a in sys.argv[1:] if not a.startswith('--')]
STEP = int(args[0]) if len(args) > 0 else 100
NSTEP = int(args[1]) if len(args) > 1 else 12
HEADED = '--headed' in sys.argv
if '--slow' in sys.argv:
    C.DRAG_STEP = 20
    C.DRAG_WAIT = int(sys.argv[sys.argv.index('--slow') + 1])
DBG = os.path.join(C.OUT, '_probe')
os.makedirs(DBG, exist_ok=True)

data = json.load(open(C.STAGES_JSON, encoding='utf-8'))
st = data['stades'][0]
t0 = time.time()
pw, ctx, page = C.start_browser(headless=not HEADED)
try:
    ok = C.wait_logged_in(page, C.fr_url(st['catalogue_voka']), timeout_s=60)
    print('logged/model page:', ok, page.url[:80], f'{time.time()-t0:.0f}s')
    box = C.viewer_box(page)
    loaded = C.wait_model_loaded(page, box, timeout_s=300)
    print('loaded:', loaded, f'{time.time()-t0:.0f}s')
    inner = [f for f in page.frames if 'unity/viewer' in f.url]
    print('inner frames:', [f.url for f in inner])
    if inner:
        fr = inner[0]
        els = fr.evaluate("""() => [...document.querySelectorAll('button, [role=button], img, svg')].map(b => {
            const r = b.getBoundingClientRect(); return {tag: b.tagName, cls: (b.className.baseVal !== undefined ? b.className.baseVal : b.className), id: b.id,
            aria: b.getAttribute('aria-label'), title: b.title, text: (b.textContent||'').trim().slice(0,30), x: Math.round(r.x), y: Math.round(r.y), w: Math.round(r.width), h: Math.round(r.height)}
        }).filter(e => e.w > 0 && e.h > 0 && e.y < 200)""")
        for e in els:
            print('  el', e)
        ast = fr.evaluate("""() => { const t = [...document.querySelectorAll('*')].filter(e => e.children.length === 0 && /Astuces/.test(e.textContent));
            return t.map(e => { const r = e.getBoundingClientRect(); return {tag: e.tagName, cls: e.className, x: r.x, y: r.y}; }); }""")
        print('  astuces:', ast)
        canv = fr.evaluate("""() => [...document.querySelectorAll('canvas')].map(c => { const r = c.getBoundingClientRect(); return {w: c.width, h: c.height, x: r.x, y: r.y, cw: r.width, ch: r.height, id: c.id}; })""")
        print('  canvas:', canv)
    page.screenshot(path=os.path.join(DBG, 'full.png'))
    cx, cy = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
    m0 = C.model_mask(C.shot(page, box)); prev = m0
    print('aire0', int(m0.sum()))
    for k in range(1, NSTEP + 1):
        C.drag(page, cx, cy, STEP, 0)
        img = C.shot(page, box); m = C.model_mask(img)
        img.resize((img.width // 3, img.height // 3)).save(os.path.join(DBG, f'p_{k * STEP:04d}.png'))
        print(f'  {k * STEP:5d} px  IoU/prev={C.similarity(prev, m):.3f}  IoU/0={C.similarity(m0, m):.3f}  aire={int(m.sum())}', flush=True)
        prev = m
    inner = [f for f in page.frames if 'unity/viewer' in f.url]
    if inner:
        els = inner[0].evaluate("""() => [...document.querySelectorAll('button, [role=button]')].map(b => {
            const r = b.getBoundingClientRect(); return {cls: b.className, aria: b.getAttribute('aria-label'), title: b.title, text: (b.textContent||'').trim().slice(0,30), x: Math.round(r.x), y: Math.round(r.y), w: Math.round(r.width), h: Math.round(r.height)}
        }).filter(e => e.w > 0 && e.h > 0)""")
        for e in els:
            print('  btn-after', e)
        page.screenshot(path=os.path.join(DBG, 'after.png'))
finally:
    ctx.close(); pw.stop()
