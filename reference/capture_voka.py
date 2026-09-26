"""
capture_voka.py — module « qui ne fait que regarder » : pilote le viewer VOKA (Unity WebGL) dans un
Chromium Playwright (headless par défaut), tourne chaque embryon J28→J91 sous des angles fixes et
enregistre les vues en PNG. Aucun maillage n'est extrait : ce sont des captures d'écran, à usage de
référence privée (forme de l'enveloppe, proportions, aspect de la peau).

Première utilisation : lancer avec --headed, se connecter soi-même à catalog.voka.io dans la fenêtre
(le profil est conservé dans .voka_profile/, la connexion n'est demandée qu'une fois), puis les
lancements suivants peuvent être headless.

Usage :
  python capture_voka.py --headed             # 1re fois (connexion), tous les stades
  python capture_voka.py --stages J28 J56     # certains stades
  python capture_voka.py --calibrate          # (re)mesure les pixels par tour (automatique la 1re fois)
  python capture_voka.py --az-step 30 --elev 0 60 -60
Sorties : captures/Jxx/Jxx_azAAA_elSEE.png + captures/Jxx/planche.png + captures/Jxx/meta.json

Structure du viewer (observée) : page Angular > iframe « unityIframe » (iframe.voka.io/model : barre
d'outils droite, boutons réinitialiser/annuler/rétablir en haut à gauche, bulle « Astuces »)
> iframe unity/viewer/index.html (canvas Unity). Un glissement lent (20 px / 30 ms) fait tourner ;
un glissement brusque = clic de sélection. Sensibilité ≈ 2250 px par tour à 60 i/s, d'où des drags segmentés.
"""
import argparse, json, os, sys, time, io
import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
STAGES_JSON = os.path.join(HERE, 'voka_stages.json')
CONFIG = os.path.join(HERE, 'capture_config.json')
PROFILE = os.path.join(HERE, '.voka_profile')
OUT = os.path.join(HERE, 'captures')
VIEWPORT = dict(width=1400, height=1000)
SETTLE = 0.4          # s après un geste, le temps que Unity redessine (pas d'inertie observée)
DRAG_STEP = 20        # px par déplacement de souris pendant un glissement
DRAG_WAIT = 30        # ms entre deux déplacements : Unity lit la souris une fois par image
FR = '/fr/'
EN = '/en/'


def fr_url(url):
    return url.replace(EN, FR)


def load_config():
    if os.path.exists(CONFIG):
        return json.load(open(CONFIG, encoding='utf-8'))
    return {}


def save_config(cfg):
    json.dump(cfg, open(CONFIG, 'w', encoding='utf-8'), indent=1)


# ----------------------------------------------------------------------------------------- navigateur
def start_browser(headless=True):
    from playwright.sync_api import sync_playwright
    pw = sync_playwright().start()
    ctx = pw.chromium.launch_persistent_context(
        PROFILE, headless=headless, viewport=VIEWPORT, locale='fr-FR',
        args=['--use-gl=angle', '--use-angle=d3d11', '--enable-gpu-rasterization', '--ignore-gpu-blocklist'])
    # d3d11 : Unity tourne à 60 i/s même en headless (SwiftShader ne fait que 10 i/s et fausse la sensibilité souris)
    page = ctx.pages[0] if ctx.pages else ctx.new_page()
    return pw, ctx, page


def wait_logged_in(page, url, timeout_s=900):
    """Va sur la page du modèle ; si le site renvoie à l'accueil (non connecté), attend que
    l'utilisateur se connecte lui-même, puis retente la page du modèle."""
    page.goto(url, wait_until='domcontentloaded')
    t0 = time.time()
    warned = False
    last_retry = 0
    while time.time() - t0 < timeout_s:
        page.wait_for_timeout(1500)
        if '/models/' in page.url and page.locator('iframe[name=unityIframe]').count() > 0:
            return True
        if not warned:
            print('\n>>> Non connecté : connectez-vous à VOKA dans la fenêtre du navigateur '
                  '(bouton « Se connecter »). Le script reprend tout seul ensuite.\n', flush=True)
            warned = True
        if '/models/' not in page.url and time.time() - last_retry > 12:
            last_retry = time.time()
            try:
                page.goto(url, wait_until='domcontentloaded')
            except Exception:
                pass
    return False


def mid_frame(page):
    for f in page.frames:
        if 'iframe.voka.io/model' in f.url:
            return f
    return None


def inner_frame(page):
    for f in page.frames:
        if 'unity/viewer' in f.url:
            return f
    return None


def viewer_box(page):
    """Rectangle de l'iframe Unity, réduit pour exclure les barres d'outils superposées."""
    b = page.locator('iframe[name=unityIframe]').first.bounding_box()
    return dict(x=b['x'] + 60, y=b['y'] + 50, width=b['width'] - 60 - 70, height=b['height'] - 50 - 60)


def dismiss_popups(page):
    """Bandeau « Annuler » de la page, bulle « Astuces » de l'iframe (croix ~245 px à droite du titre)."""
    try:
        loc = page.locator('button:has-text("Annuler")')
        if loc.count() and loc.first.is_visible():
            loc.first.click(timeout=800)
    except Exception:
        pass
    for f in [mid_frame(page), inner_frame(page), page]:
        if f is None:
            continue
        try:
            t = f.get_by_text('Astuces', exact=True)
            if t.count() and t.first.is_visible():
                b = t.first.bounding_box()
                page.mouse.click(b['x'] + 245, b['y'] + b['height'] / 2)
                page.wait_for_timeout(300)
                return
        except Exception:
            pass


def collapse_thumbnails(page):
    """Replie la colonne des vignettes de stades (bouton « Modèles » de la barre gauche) si elle est ouverte."""
    try:
        if page.get_by_text('Jour 35', exact=True).first.is_visible():
            page.locator('button:has-text("Modèles")').first.click(timeout=800)
            page.wait_for_timeout(500)
    except Exception:
        pass


def wait_model_loaded(page, box, timeout_s=300):
    """Attend la fin du chargement Unity : canvas dimensionné, écran de progression disparu,
    puis silhouette non vide et stable sur deux prises consécutives."""
    t0 = time.time()
    prev = None
    while time.time() - t0 < timeout_s:
        fr = inner_frame(page)
        ready = False
        if fr:
            try:
                ready = fr.evaluate("""() => { const c = document.querySelector('canvas');
                    const txt = document.body.innerText || '';
                    return !!c && c.width > 300 && !/%|Finalisation|Chargement|Loading/i.test(txt); }""")
            except Exception:
                ready = False
        if ready:
            m = model_mask(shot(page, box))
            area = int(m.sum())
            if prev is not None and area > 0.01 * m.size and abs(area - prev) < 0.02 * area:
                page.wait_for_timeout(500)
                return True
            prev = area
        page.wait_for_timeout(1000)
    print('!! modèle non chargé après %ds' % timeout_s)
    return False


def reset_view(page):
    """Bouton « réinitialiser la vue » = premier .button (haut-gauche) de l'iframe intermédiaire."""
    f = mid_frame(page)
    try:
        if f is not None:
            f.locator('.button').first.click(timeout=1000)
            page.wait_for_timeout(700)
            return
    except Exception:
        pass
    b = page.locator('iframe[name=unityIframe]').first.bounding_box()
    page.mouse.click(b['x'] + 82, b['y'] + 28)
    page.wait_for_timeout(700)


def drag(page, cx, cy, dx, dy):
    """Glissement progressif (rotation orbitale). Le curseur est posé au point de départ 150 ms AVANT
    l'appui : sinon Unity compte le saut du curseur dans le geste (chaque drag ajoutait une rotation parasite)."""
    page.mouse.move(cx, cy)
    page.wait_for_timeout(150)
    page.mouse.down()
    page.wait_for_timeout(40)
    n = max(1, int(max(abs(dx), abs(dy)) / DRAG_STEP))
    for i in range(1, n + 1):
        page.mouse.move(cx + dx * i / n, cy + dy * i / n)
        page.wait_for_timeout(DRAG_WAIT)
    page.wait_for_timeout(40)
    page.mouse.up()
    page.wait_for_timeout(int(SETTLE * 1000))


def drag_total(page, box, dx=0, dy=0):
    """Glissement de grande amplitude découpé en segments qui restent dans la zone 3D."""
    cx, cy = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
    segx = box['width'] * 0.6
    segy = box['height'] * 0.6
    while abs(dx) > 1 or abs(dy) > 1:
        sx = max(-segx, min(segx, dx))
        sy = max(-segy, min(segy, dy))
        drag(page, cx - sx / 2, cy - sy / 2, sx, sy)
        dx -= sx
        dy -= sy


def shot(page, box):
    png = page.screenshot(clip=box)
    return Image.open(io.BytesIO(png)).convert('RGB')


def model_mask(img, tol=18):
    """Masque du modèle = pixels différents de la couleur de fond (médiane des bords)."""
    a = np.asarray(img).astype(np.int16)
    border = np.concatenate([a[0], a[-1], a[:, 0], a[:, -1]])
    bg = np.median(border, axis=0)
    m = (np.abs(a - bg).max(axis=2) > tol)
    m[-110:, -380:] = False    # coin bas-droit : bulle « Astuces » éventuelle
    return m


def similarity(m1, m2):
    inter = np.logical_and(m1, m2).sum()
    uni = np.logical_or(m1, m2).sum()
    return inter / max(uni, 1)


def prepare(page, url):
    """Ouvre la page du modèle, attend le chargement, range l'interface. Renvoie la zone 3D."""
    if not wait_logged_in(page, fr_url(url)):
        return None
    page.wait_for_timeout(2000)
    box = viewer_box(page)
    if not wait_model_loaded(page, box):
        return None
    dismiss_popups(page)
    collapse_thumbnails(page)
    return box


# ----------------------------------------------------------------------------------------- calibrage
def calibrate(page, box, step_px=100, max_px=4000):
    """Tourne par pas cumulés et cherche le glissement qui ramène la silhouette initiale : px par 360°."""
    reset_view(page)
    dbg = os.path.join(OUT, '_calib')
    os.makedirs(dbg, exist_ok=True)
    page.screenshot(path=os.path.join(dbg, 'page_full.png'))
    img0 = shot(page, box)
    img0.save(os.path.join(dbg, 'calib_0000.png'))
    ref = model_mask(img0)
    sims, areas = [], []
    for k in range(1, max_px // step_px + 1):
        drag_total(page, box, dx=step_px)
        img = shot(page, box)
        if k % 5 == 0:
            img.resize((img.width // 2, img.height // 2)).save(os.path.join(dbg, f'calib_{k * step_px:04d}.png'))
        m = model_mask(img)
        sims.append(similarity(ref, m))
        areas.append(int(m.sum()))
        print(f'  calib {k * step_px:5d} px  IoU={sims[-1]:.3f}  aire={areas[-1]}', flush=True)
        # arrêt : revenu à la silhouette initiale après s'en être éloigné, et déjà passé le maximum
        if k > 10 and min(sims) < 0.5 and max(sims[-6:]) > 0.85 and sims[-1] < max(sims[-6:]) - 0.05:
            break
    sims = np.array(sims)
    started = int(np.argmax(sims < 0.5))
    k = started + int(np.argmax(sims[started:]))
    px360 = (k + 1) * step_px
    print(f'>>> ~{px360} px par tour (IoU {sims[k]:.3f})')
    return px360, sims.tolist(), areas


# ----------------------------------------------------------------------------------------- captures
def capture_stage(page, st, cfg, az_step, elevs, force=False):
    sid = st['id']
    d = os.path.join(OUT, sid)
    os.makedirs(d, exist_ok=True)
    print(f'=== {sid} ({st["carnegie"]}) {st["catalogue_voka"]}')
    box = prepare(page, st['catalogue_voka'])
    if box is None:
        print('!! page ou modèle indisponible, stade sauté')
        return
    px360 = cfg['px_per_360']
    px_el = cfg.get('px_per_90_elev', px360 / 4)
    views = []
    for el in elevs:
        for az in range(0, 360, az_step):
            fn = f'{sid}_az{az:03d}_el{el:+03d}.png'
            fp = os.path.join(d, fn)
            if os.path.exists(fp) and not force:
                views.append(dict(file=fn, az=az, el=el))
                continue
            reset_view(page)
            dismiss_popups(page)
            if el:
                drag_total(page, box, dy=-el / 90 * px_el)   # glisser vers le haut = voir par-dessus
            if az:
                drag_total(page, box, dx=az / 360 * px360)
            shot(page, box).save(fp)
            views.append(dict(file=fn, az=az, el=el))
            print(f'  {fn}', flush=True)
    meta = dict(stage=st, views=views, box=box, px_per_360=px360, viewport=VIEWPORT,
                note='az = rotation horizontale depuis la vue par défaut (0 = vue latérale VOKA) ; '
                     'el > 0 = vu de dessus. Fond uniforme -> silhouettes par silhouettes.py')
    json.dump(meta, open(os.path.join(d, 'meta.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    planche(d, views, sid)


def planche(d, views, sid, thumb=260):
    ims = []
    for v in views:
        im = Image.open(os.path.join(d, v['file'])).convert('RGB')
        im.thumbnail((thumb, thumb))
        ims.append((im, v))
    if not ims:
        return
    els = sorted(set(v['el'] for _, v in ims))
    azs = sorted(set(v['az'] for _, v in ims))
    sheet = Image.new('RGB', (thumb * len(azs), thumb * len(els)), (200, 200, 200))
    dr = ImageDraw.Draw(sheet)
    for im, v in ims:
        x = azs.index(v['az']) * thumb
        y = els.index(v['el']) * thumb
        sheet.paste(im, (x + (thumb - im.width) // 2, y + (thumb - im.height) // 2))
        dr.text((x + 4, y + 4), f"{sid} az{v['az']} el{v['el']}", fill=(40, 40, 40))
    sheet.save(os.path.join(d, 'planche.png'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', nargs='*', help='J28 J35 ... (défaut : tous)')
    ap.add_argument('--az-step', type=int, default=30)
    ap.add_argument('--elev', type=int, nargs='*', default=[0, 60, -60])
    ap.add_argument('--calibrate', action='store_true')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--headed', action='store_true', help='fenêtre visible (nécessaire pour la 1re connexion)')
    a = ap.parse_args()
    data = json.load(open(STAGES_JSON, encoding='utf-8'))
    stages = [s for s in data['stades'] if not a.stages or s['id'] in a.stages]
    cfg = load_config()
    pw, ctx, page = start_browser(headless=not a.headed)
    try:
        if a.calibrate or 'px_per_360' not in cfg:
            st = stages[0]
            box = prepare(page, st['catalogue_voka'])
            if box is None:
                sys.exit('page ou modèle indisponible')
            px360, sims, areas = calibrate(page, box)
            cfg.update(px_per_360=px360, calib_iou=sims, calib_area=areas, calib_stage=st['id'])
            save_config(cfg)
        for st in stages:
            capture_stage(page, st, cfg, a.az_step, a.elev, a.force)
    finally:
        ctx.close()
        pw.stop()


if __name__ == '__main__':
    main()
