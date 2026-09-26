# -*- coding: utf-8 -*-
"""Viewer épuré « embryon humain en croissance, par système » → embryons_3D/site_3dh/index.html (harcelon.fr/3dh).

    python embryo3d/viewer_3dh.py                 # écrit site_3dh/index.html (GLB : ../3dht/glb/, version en ligne)
    python embryo3d/viewer_3dh.py --local         # écrit site_3dh/index_local.html (GLB : ../site_3dht/glb/, pour tester ici)

Lit les manifestes de stade (<stade>/out/manifest.json) pour connaître les structures présentes ; les GLB sont ceux du
paquet 3dht (mêmes fichiers, pas de nouvel envoi). Un système = une couleur ; les brouillons du pipeline remplacés par une
reconstruction ne sont jamais affichés.
"""
import json
import os
import sys
from datetime import datetime

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)
STADES = [("CS13", "CS13.f4v", 28, 32), ("CS14", "CS14_f4v", 31, 35), ("CS15", "CS15_f4v", 35, 38), ("CS16", "CS16_f4v", 37, 42),
          ("CS17", "CS17_f4v", 42, 44), ("CS19", "CS19_f4v", 48, 51), ("CS20", "CS20_F4V", 51, 53)]
LOCAL = "--local" in sys.argv
EMBED = "--embed" in sys.argv     # viewer.html dans le paquet 3dht (GLB et topo du même dossier)
GLB_BASE = "glb/" if EMBED else ("../site_3dht/glb/" if LOCAL else "../3dht/glb/")
TOPO_BASE = "topo/" if EMBED else ("../site_3dht/topo/" if LOCAL else "../3dht/topo/")
GL_MM = {"CS13": 4.5, "CS14": 6.0, "CS15": 7.5, "CS16": 9.0, "CS17": 11.0, "CS19": 17.0, "CS20": 20.0}


def boite_glb(chemin, gl):
    """Boîte englobante (mm) du modèle complet, lue dans les accesseurs du GLB : taille connue de l'embryon,
    indépendante des structures affichées. Repli sur les proportions moyennes (0,42 × 0,69 × 1 GL) si le fichier manque."""
    try:
        import struct
        with open(chemin, "rb") as f:
            f.read(12); cl, _ = struct.unpack("<II", f.read(8)); g = json.loads(f.read(cl))
        mn, mx = [1e9] * 3, [-1e9] * 3
        for m in g["meshes"]:
            for pr in m["primitives"]:
                a = g["accessors"][pr["attributes"]["POSITION"]]
                for i in range(3):
                    mn[i] = min(mn[i], a["min"][i]); mx[i] = max(mx[i], a["max"][i])
        return {"min": [round(x, 3) for x in mn], "max": [round(x, 3) for x in mx]}
    except Exception:
        return {"min": [-0.21 * gl, -0.35 * gl, -0.5 * gl], "max": [0.21 * gl, 0.35 * gl, 0.5 * gl]}

# système → (couleur, structures du GLB, opacité)
SYSTEMES = [
    ("nerveux", "Système nerveux (segmentation fine)", "#5b8be0", {"prosencephale": "#3f6fc4", "mesencephale": "#5b8be0", "rhombencephale": "#7fa8f0", "ventricule_prosencephale": "#bfe3ff", "ventricule_mesencephale": "#bfe3ff", "ventricule_rhombencephale": "#bfe3ff", "canal_central": "#dff2ff", "moelle": "#2f5fb0", "moelle_rachidienne": "#2f5fb0", "meninges_mesenchyme_cranien": "#a9b8c8"}),
    ("snc", "Système nerveux (segmentation automatique)", "#4f7fc9", {"snc": "#4f7fc9", "ventricules": "#8fd3ff", "ganglions": "#9fc5e8"}),
    ("sens", "Yeux et otocystes", "#6aa7e8", {"yeux": "#2f4f9f", "cristallins": "#dfe9ff", "vesicules_otiques": "#6aa7e8"}),
    ("coeur", "Cœur", "#b0202a", {"coeur_detoure": "#b0202a", "myocarde": "#b0202a", "cavites_cardiaques": "#f08a8a"}),
    ("arteres", "Artères", "#d62828", {"vaisseaux_aorte": "#d62828"}),
    ("veines", "Veines", "#2f6fd6", {"vaisseaux_cardinales": "#2f6fd6", "vaisseaux_ombilicaux": "#3d8bff", "vaisseaux_vitellins": "#1f4fd6", "vaisseaux_veines": "#2f6fd6"}),
    ("digestif", "Tube digestif", "#f2c14e", {"digestif_oesophage": "#f2c14e", "digestif_estomac": "#f2c14e", "digestif_duodenum": "#f2c14e", "digestif_intestin_moyen": "#f2c14e", "digestif_intestin_posterieur": "#f2c14e"}),
    ("foie", "Foie", "#8c5a2b", {"foie": "#8c5a2b"}),
    ("colonne", "Squelette axial : cartilages, arcs, côtes, notochorde, axe", "#f1ead6", {"axe_vertebral": "#ff9f43", "etages_vertebraux": "#f1ead6", "somites_video": "#ffd27f", "etages_manuels": "#7fe07f", "etages_manuels_dorsal": "#4fc94f", "reper1": "#bfffbf", "reper2": "#bfffbf", "notochorde": "#ffffff", "squelette_axial_cartilage": "#f1ead6", "chondrocrane": "#e8e0c8", "cartilage_autre": "#d9d0b8", "corps_vertebraux": "#f5f0e0", "arcs_neuraux": "#e6dcc0", "cotes": "#d8ccb0"}),
    ("muscles", "Somites", "#9b59b6", {"somites": "#9b59b6"}),
    ("appendiculaire", "Appendiculaire (membres)", "#e0a890", {"membres": "#e0a890"}),
    ("derme", "Derme (peau, cordon)", "#e8b8a0", {"enveloppe": "#e8b8a0", "cordon_ombilical": "#d8a088"}),
]
TRANSLUCIDES = {"enveloppe": 0.16}

stades = []
for cs, d, j0, j1 in STADES:
    m = json.load(open(os.path.join(d, "out", "manifest.json"), encoding="utf-8"))
    noms = [x["name"] for x in m["structures"] if x.get("confiance") != "faible"]   # pièces réellement affichables (GLB allégé)
    musc = []   # superpositions : topographie (axe vertébral, étages) seulement si fiable
    mt = os.path.join(d, "out", "topographie", "manifest_topographie.json")
    if os.path.exists(mt):
        t = json.load(open(mt, encoding="utf-8")); ct = t.get("controle") or {}
        fiable = (ct.get("etages") or {}).get("fiabilite", {}).get("fiable") and not ct.get("avertissement")
        for x in t["segments"]:
            if fiable or x.get("confiance") == "utilisateur":     # marques manuelles de l'utilisateur : toujours affichées
                musc.append({"nom": x["name"], "url": TOPO_BASE + os.path.basename(x["file"]), "confiance": "bonne"})   # paquet site : fichiers aplatis dans topo/
    glb = GLB_BASE + f"{cs}_all.glb"   # modèle complet : tout le travail du pipeline
    gl = GL_MM.get(cs, 20.0)
    stades.append({"id": cs, "j0": j0, "j1": j1, "mid": (j0 + j1) / 2, "glb": glb, "structures": noms, "musc": musc, "gl_mm": gl,
                   "box": boite_glb(os.path.join("embryons_3D", "site_3dht", "glb", f"{cs}_all.glb"), gl)})

CARTILAGE_OK = ["CS17", "CS20"]   # arcs et côtes déclarés bonne / moyenne par la session base ; ailleurs = fragments automatiques
VERSION_FILE = "viewer_version.txt" if EMBED else "version.txt"
data = {"stades": stades, "systemes": [{"id": i, "lbl": l, "col": c, "structs": st} for i, l, c, st in SYSTEMES], "cartilage_ok": CARTILAGE_OK, "version_file": VERSION_FILE,
        "translucides": TRANSLUCIDES, "version": datetime.now().strftime("%Y%m%d-%H%M%S")}

HTML = r"""<!DOCTYPE html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="robots" content="noindex">
<title>Embryon humain · croissance par système</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
:root{--bg:#0e1013;--panel:#171a1f;--line:#262a31;--ink:#eceae4;--mut:#8f948f;--acc:#ff9f43}
*{box-sizing:border-box}html,body{height:100%;margin:0;background:var(--bg);color:var(--ink);font:14px/1.4 Inter,system-ui,sans-serif}
#app{display:grid;grid-template-columns:260px 1fr;grid-template-rows:1fr auto;height:100%}
aside{grid-row:1/3;background:var(--panel);border-right:1px solid var(--line);padding:18px 16px;display:flex;flex-direction:column;gap:14px;overflow:auto}
h1{font-size:15px;margin:0;font-weight:700;letter-spacing:-.01em}h1 small{display:block;color:var(--mut);font-weight:400;font-size:12px;margin-top:2px}
.stage{font-size:30px;font-weight:700;line-height:1;margin-top:4px}.stage small{display:block;font-size:12px;color:var(--mut);font-weight:500;margin-top:4px}
.sys{display:flex;flex-direction:column;gap:6px}.sys button{display:flex;align-items:center;gap:10px;background:transparent;border:1px solid var(--line);color:var(--ink);padding:8px 10px;border-radius:9px;cursor:pointer;text-align:left;font:inherit}
.sys button i{width:12px;height:12px;border-radius:50%;flex:0 0 auto}.sys button.off{opacity:.35}.sys button .n{margin-left:auto;color:var(--mut);font-size:11px}
.row{display:flex;gap:6px;flex-wrap:wrap}.row button,.bar button{background:transparent;border:1px solid var(--line);color:var(--ink);padding:6px 10px;border-radius:8px;cursor:pointer;font:inherit;font-size:13px}
.bar button.on,.row button.on{background:var(--acc);border-color:var(--acc);color:#111}
label.chk{display:flex;align-items:center;gap:6px;font-size:13px;color:var(--mut)}
main{position:relative;min-height:0}canvas{width:100%;height:100%;display:block}
#msg{position:absolute;left:14px;top:12px;font-size:12px;color:var(--mut);background:rgba(0,0,0,.4);padding:4px 8px;border-radius:6px;pointer-events:none}
.bar{grid-column:2;display:flex;align-items:center;gap:12px;padding:10px 16px;border-top:1px solid var(--line);background:var(--panel)}
.tl{flex:1;position:relative;height:38px}.tl input{position:absolute;left:0;top:14px;width:100%;margin:0;-webkit-appearance:none;background:transparent}
.tl input::-webkit-slider-runnable-track{height:4px;background:var(--line);border-radius:2px}.tl input::-webkit-slider-thumb{-webkit-appearance:none;width:16px;height:16px;border-radius:50%;background:var(--acc);margin-top:-6px;cursor:grab}
.tl .tk{position:absolute;top:0;transform:translateX(-50%);font-size:11px;color:var(--mut)}.tl .tk::after{content:'';display:block;width:1px;height:6px;background:var(--mut);margin:1px auto 0}
.bar .day{font-weight:600;min-width:64px;text-align:right}
#tab,#fs,#badge{display:none}#app{position:relative}
#badge{position:absolute;right:10px;top:8px;font-weight:700;font-size:18px;background:rgba(0,0,0,.45);padding:3px 9px;border-radius:8px;pointer-events:none}
@media (max-width:800px),(max-height:520px){#app{grid-template-columns:1fr;grid-template-rows:1fr auto}main{grid-row:1;grid-column:1}.bar{grid-column:1;grid-row:2;padding:6px 8px;gap:8px}
 aside{position:absolute;left:0;top:0;bottom:0;width:min(78vw,300px);z-index:5;grid-row:auto;padding:12px;gap:10px;transform:translateX(-100%);transition:transform .22s;box-shadow:0 0 30px rgba(0,0,0,.6)}
 body.drawer aside{transform:none}h1 small{display:none}.stage{font-size:22px}.sys button{padding:6px 9px;font-size:13px}
 #tab,#fs{display:inline-flex;align-items:center;justify-content:center;min-width:38px}#badge{display:block}.bar .day{min-width:44px}.tl .tk{display:none}.tl .tk:nth-child(odd){display:block}}
</style></head><body><div id="app">
<aside>
 <h1>Embryon humain <small>croissance par système · CS13 → CS20 · première version (24/09)</small></h1>
 <div class="stage"><span id="st">—</span><small id="stInfo">chargement…</small></div>
 <div class="sys" id="sys"></div>
 <div class="row"><button id="all">tout</button><button id="none">rien</button><button id="fit">recadrer</button></div>
 <label class="chk" id="morphLbl" style="display:none"><input type="checkbox" id="morph" checked> morphing continu (clés de la scène Blender)</label>
 <label class="chk"><input type="checkbox" id="autozoom"> zoom interpolé sur la croissance (décoché : échelle constante, caméra fixe)</label>
 <label class="chk" id="tubesLbl"><input type="checkbox" id="tubes" checked> tubes interpolés entre stades (digestif, aortes)</label>
 <label class="chk"><input type="checkbox" id="loop" checked> lecture en boucle</label>
 <label class="chk">vitesse <select id="speed"><option value="1">1 j/s</option><option value="2" selected>2 j/s</option><option value="4">4 j/s</option></select></label>
 <div style="margin-top:auto;font-size:11px;color:var(--mut)">glisser = orbite · molette = zoom · clic droit = déplacer<br>version <span id="ver"></span></div>
</aside>
<main><canvas id="gl"></canvas><div id="msg"></div><div id="badge">—</div></main>
<div class="bar"><button id="tab" title="calques">☰</button><button id="play" class="on">❚❚</button><div class="tl" id="tl"><input type="range" id="day" min="26" max="55" step="0.05" value="28"></div><div class="day" id="dayLbl">J28</div><button id="fs" title="plein écran">⛶</button></div>
</div>
<script id="data" type="application/json">__DATA__</script>
<script type="importmap">{"imports":{"three":"https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/","three-mesh-bvh":"https://cdn.jsdelivr.net/npm/three-mesh-bvh@0.7.8/build/index.module.js"}}</script>
<script type="module">
import * as THREE from 'three';import {OrbitControls} from 'three/addons/controls/OrbitControls.js';import {GLTFLoader} from 'three/addons/loaders/GLTFLoader.js';import {PLYLoader} from 'three/addons/loaders/PLYLoader.js';
let MeshBVH=null;try{({MeshBVH}=await import('three-mesh-bvh'));}catch(e){}
const D=JSON.parse(document.getElementById('data').textContent);const $=id=>document.getElementById(id);$('ver').textContent=D.version;
setInterval(async()=>{try{const r=await fetch(D.version_file+'?t='+Date.now(),{cache:'no-store'});if(r.ok){const v=(await r.text()).trim();if(v&&v!==D.version)location.reload();}}catch(e){}},20000);
const ST=D.stades;const SYS=D.systemes;const colorOf={},sysOf={};SYS.forEach(s=>Object.entries(s.structs).forEach(([n,c])=>{colorOf[n]=c;sysOf[n]=s.id;}));

// lecture continue (consigne utilisateur 24/09 15:45 : « les choses qui apparaissent d'un coup, on les enlève ») :
// pendant la lecture, seules les pièces présentes aux 7 stades sont affichées ; en pause, tout le stade.
const COMMUN=new Set([...ST[0].structures.filter(n=>ST.every(s=>s.structures.includes(n))),...ST[0].musc.map(m=>m.nom).filter(n=>ST.every(s=>s.musc.some(m=>m.nom===n)))]);
const SCRATCH=new Set(['coeur','tube_digestif','vaisseaux','cavite_pericardique','muscles_membres']); // snc : interdit tant qu'il n'est pas segmenté (consigne 24/09)
// ---- scène
const canvas=$('gl');const R=new THREE.WebGLRenderer({canvas,antialias:true,preserveDrawingBuffer:true});R.setPixelRatio(Math.min(devicePixelRatio,2));
const scene=new THREE.Scene();scene.background=new THREE.Color(0x0e1013);
scene.add(new THREE.HemisphereLight(0xffffff,0x334455,2.4));const key=new THREE.DirectionalLight(0xffffff,2.6);key.position.set(-3,2,4);scene.add(key);const fill=new THREE.DirectionalLight(0xffffff,1.0);fill.position.set(3,-2,-2);scene.add(fill);
const cam=new THREE.PerspectiveCamera(32,1,0.05,500);{const t=THREE.MathUtils.degToRad(15);const v=new THREE.Vector3(Math.cos(t),-Math.sin(t),0).normalize();const Dd=new THREE.Vector3(0,1,0);const Dp=Dd.clone().sub(v.clone().multiplyScalar(Dd.dot(v))).normalize();const w=new THREE.Vector3().crossVectors(v,Dp);const sd=new THREE.Vector2(-0.6,0.8).normalize();cam.up.copy(Dp.clone().multiplyScalar(sd.y).add(w.clone().multiplyScalar(-sd.x)).normalize());}
const ctl=new OrbitControls(cam,canvas);ctl.enableDamping=true;ctl.dampingFactor=0.08;ctl.rotateSpeed=0.75;ctl.zoomSpeed=0.8;ctl.panSpeed=0.7;ctl.screenSpacePanning=true;
let lastAspect=0;function resize(){const w=canvas.clientWidth,h=canvas.clientHeight;if(!w||!h)return;R.setSize(w,h,false);cam.aspect=w/h;cam.updateProjectionMatrix();
 const big=Math.abs(cam.aspect-lastAspect)>0.01*Math.max(1,lastAspect);lastAspect=cam.aspect;if(big&&typeof fit==='function'&&fit.ready&&!userMoved)fit(false);}new ResizeObserver(resize).observe(canvas);resize();
// vue de référence : côté gauche (vers +X), dos incliné de 15° vers la caméra, dos en haut à gauche (repère X droite, Y dorsal, Z haut)
function refPose(){const t=THREE.MathUtils.degToRad(15);const v=new THREE.Vector3(Math.cos(t),-Math.sin(t),0).normalize();const Dd=new THREE.Vector3(0,1,0);
 const Dp=Dd.clone().sub(v.clone().multiplyScalar(Dd.dot(v))).normalize();const w=new THREE.Vector3().crossVectors(v,Dp);const sd=new THREE.Vector2(-0.6,0.8).normalize();
 const up=Dp.clone().multiplyScalar(sd.y).add(w.clone().multiplyScalar(-sd.x)).normalize();return {v,up};}
let userMoved=false;ctl.addEventListener('start',()=>{userMoved=true;camGoal=null;});
// ---- état
const OFF_DEFAUT=['snc','muscles','appendiculaire','derme'];const on={};SYS.forEach(s=>on[s.id]=!OFF_DEFAUT.includes(s.id));
const loaded={},loading={},plyCache={};let root=null,cur=null,muscGroup=null;
let day=28,playing=true,last=performance.now();
const stageAt=d=>{let r=ST[0];for(const s of ST)if(s.mid<=d)r=s;return r;};const nextOf=s=>{const i=ST.indexOf(s);return i<ST.length-1?ST[i+1]:null;};
function load(s,cb){if(loaded[s.id]){cb&&cb();return;}if(loading[s.id]){return;}loading[s.id]=true;$('msg').textContent='chargement '+s.id+'…';
 new GLTFLoader().load(s.glb,g=>{g.scene.traverse(o=>{if(o.isMesh){if(!o.geometry.attributes.normal)o.geometry.computeVertexNormals();o.material=new THREE.MeshStandardMaterial({color:new THREE.Color(colorOf[o.name]||'#777'),roughness:0.55,metalness:0,side:THREE.DoubleSide});
   if(D.translucides[o.name]!=null){o.material.transparent=true;o.material.opacity=D.translucides[o.name];o.material.depthWrite=false;}}});loaded[s.id]=g.scene;loading[s.id]=false;$('msg').textContent='';cb&&cb();},
  x=>{if(x.total)$('msg').textContent=`chargement ${s.id} : ${Math.round(100*x.loaded/x.total)} %`;},()=>{loading[s.id]=false;$('msg').textContent='échec '+s.id;});
 s.musc.forEach(m=>{if(!plyCache[m.url]){plyCache[m.url]='…';new PLYLoader().load(m.url,geo=>{geo.computeVertexNormals();plyCache[m.url]=geo;});}});}
function show(s){if(cur===s.id)return;if(!loaded[s.id]){load(s);return;}cur=s.id;if(root)scene.remove(root);root=loaded[s.id];scene.add(root);
 if(muscGroup)scene.remove(muscGroup);muscGroup=new THREE.Group();s.musc.forEach(m=>{const g=plyCache[m.url];if(g&&g.attributes){const mesh=new THREE.Mesh(g,new THREE.MeshStandardMaterial({color:new THREE.Color(colorOf[m.nom]||'#ff9f43'),roughness:0.55}));mesh.name=m.nom;muscGroup.add(mesh);}});scene.add(muscGroup);
 applyVis();$('st').textContent=s.id;$('badge').textContent=s.id;const n=nextOf(s);if(n)load(n);if(!fitted){fit();fitted=true;}updateInfo();}
let fitted=false,camGoal=null;const camVel={p:new THREE.Vector3(),t:new THREE.Vector3()};
const TUBE_SYS=n=>/^aorte/.test(n)?'arteres':'digestif';const TUBE_COL=n=>/^aorte/.test(n)?'#d62828':'#f2c14e';
function tubeMasque(name,s){if(!D.tubes||!$('tubes').checked||!s)return false;const P=D.tubes.presence;
 if(name==='vaisseaux_aorte')return D.tubes.noms.some(n=>/^aorte/.test(n)&&P[n]&&P[n][s.id]);
 const m=/^digestif_(.+)$/.exec(name);return !!(m&&P[m[1]]&&P[m[1]][s.id]);}
function visPred(rt){let fin=false,recr=false;const sCur=stageAt(day);const continu=playing;if(rt)rt.traverse(o=>{if(o.isMesh&&/^(prosencephale|rhombencephale|moelle|moelle_rachidienne)$/.test(o.name))fin=true;if(o.isMesh&&/^(corps_vertebraux|arcs_neuraux)$/.test(o.name))recr=true;});return o=>{if(!o.isMesh)return false;if(continu&&!COMMUN.has(o.name)&&!/^tube:/.test(o.name))return false;const sys=sysOf[o.name];if(fin&&['snc','ventricules','ganglions'].includes(o.name))return false;if(recr&&o.name==='squelette_axial_cartilage')return false;if(tubeMasque(o.name,sCur))return false;return !SCRATCH.has(o.name)&&!!sys&&on[sys]!==false;};}
function applyVis(){const pr=visPred(root);const f=o=>{if(o.isMesh)o.visible=pr(o);};if(root)root.traverse(f);if(muscGroup)muscGroup.traverse(f);}
// boîte des structures visibles (enveloppe exclue) ; fonctionne aussi pour un stade chargé mais pas encore affiché
function boxOf(rt,mg){const b=new THREE.Box3(),t=new THREE.Box3();let any=false;if(!rt)return b;rt.updateMatrixWorld(true);const pr=visPred(rt);const add=o=>{if(!o.geometry.boundingBox)o.geometry.computeBoundingBox();t.copy(o.geometry.boundingBox).applyMatrix4(o.matrixWorld);b.union(t);any=true;};
 rt.traverse(o=>{if(o.isMesh&&pr(o)&&o.name!=='enveloppe')add(o);});if(mg)mg.traverse(o=>{if(o.isMesh&&o.visible)add(o);});if(!any)rt.traverse(o=>{if(o.isMesh)add(o);});return b;}
function visibleBox(){return boxOf(root,muscGroup);}
// cadrage sur la taille connue de l'embryon : boîte du modèle complet embarquée dans les données (jamais les structures visibles,
// dont les toggles ne bougent donc plus la caméra). Zoom auto = boîte du stade ; sinon = boîte commune à tous les stades (échelle constante).
const UNION=new THREE.Box3();ST.forEach(s=>UNION.union(new THREE.Box3(new THREE.Vector3(...s.box.min),new THREE.Vector3(...s.box.max))));
function boxAt(d){let a=ST[0],b=ST[0];for(const s of ST){if(s.mid<=d){a=s;}else{b=s;break;}}if(b===a||d<=ST[0].mid){const x=d<=ST[0].mid?ST[0]:a;return new THREE.Box3(new THREE.Vector3(...x.box.min),new THREE.Vector3(...x.box.max));}
 const k=Math.max(0,Math.min(1,(d-a.mid)/(b.mid-a.mid)));const mn=new THREE.Vector3(...a.box.min).lerp(new THREE.Vector3(...b.box.min),k),mx=new THREE.Vector3(...a.box.max).lerp(new THREE.Vector3(...b.box.max),k);return new THREE.Box3(mn,mx);}
function knownBox(s){return $('autozoom').checked?boxAt(day):UNION.clone();}
function framing(b,keepDir){const c=b.getCenter(new THREE.Vector3()),sz=b.getSize(new THREE.Vector3());const r=sz.length()/2+0.15;
 const fv=THREE.MathUtils.degToRad(cam.fov)/2,fh=Math.atan(Math.tan(fv)*Math.max(0.3,Math.min(5,cam.aspect||1)));
 let dir;if(keepDir){dir=cam.position.clone().sub(ctl.target);if(dir.lengthSq()<1e-6)dir=refPose().v.clone().negate();dir.normalize();}else dir=refPose().v.clone().negate();
 // ajustement serré : les 8 coins de la boîte projetés dans le champ (largeur ET hauteur), marge 6 %
 const fwd=dir.clone().negate(),right=new THREE.Vector3().crossVectors(fwd,cam.up).normalize(),upv=new THREE.Vector3().crossVectors(right,fwd).normalize();const tv=Math.tan(fv),th=Math.tan(fh);let need=0;
 for(let i=0;i<8;i++){const q=new THREE.Vector3(i&1?b.max.x:b.min.x,i&2?b.max.y:b.min.y,i&4?b.max.z:b.min.z).sub(c);const z=q.dot(dir);need=Math.max(need,Math.abs(q.dot(right))/th+z,Math.abs(q.dot(upv))/tv+z);}
 const dist=Math.max(need*1.06+0.1,r/Math.sin(Math.min(fv,fh))*0.55);
 return {pos:c.clone().add(dir.multiplyScalar(dist)),tgt:c,dist};}
function fit(doux){fit.ready=true;if(!root)return;const s=stageAt(day);const g=framing(knownBox(s),doux&&userMoved);camGoal=g;cam.near=Math.max(0.01,g.dist/300);cam.far=Math.max(cam.far,g.dist*6);cam.updateProjectionMatrix();
 if(!doux){cam.position.copy(g.pos);ctl.target.copy(g.tgt);camGoal=null;camVel.p.set(0,0,0);camVel.t.set(0,0,0);ctl.update();}}
// glissement critique (SmoothDamp) : vitesse continue, cible re-visable à tout moment sans à-coup
function sdamp(cur,tgt,vel,st,dt){const o=2/st,x=o*dt,e=1/(1+x+0.48*x*x+0.235*x*x*x);const ch=cur.clone().sub(tgt);const tmp=vel.clone().addScaledVector(ch,o).multiplyScalar(dt);vel.sub(tmp.clone().multiplyScalar(o)).multiplyScalar(e);cur.copy(tgt).add(ch.add(tmp).multiplyScalar(e));}
$('fit').onclick=()=>{userMoved=false;fit(true);};$('autozoom').onchange=()=>{fit(true);};
// ---- systèmes
$('sys').innerHTML=SYS.map(s=>`<button data-id="${s.id}" class="${on[s.id]?'':'off'}"><i style="background:${s.col}"></i>${s.lbl}<span class="n" data-n="${s.id}"></span></button>`).join('');
document.querySelectorAll('#sys button').forEach(b=>b.onclick=()=>{on[b.dataset.id]=!on[b.dataset.id];b.classList.toggle('off',!on[b.dataset.id]);applyVis();updateMorphWeb(true);});
$('all').onclick=()=>{SYS.forEach(s=>on[s.id]=true);document.querySelectorAll('#sys button').forEach(b=>b.classList.remove('off'));applyVis();updateMorphWeb(true);};
$('none').onclick=()=>{SYS.forEach(s=>on[s.id]=false);document.querySelectorAll('#sys button').forEach(b=>b.classList.add('off'));applyVis();updateMorphWeb(true);};
function updateInfo(){const s=stageAt(day);if(!root){return;}SYS.forEach(sy=>{let n=0;root.traverse(o=>{if(o.isMesh&&sysOf[o.name]===sy.id&&!SCRATCH.has(o.name))n++;});if(sy.id==='muscles'&&muscGroup)n+=muscGroup.children.length;const el=document.querySelector(`[data-n="${sy.id}"]`);if(el)el.textContent=n?n:'—';});
 const nx=nextOf(s);$('stInfo').textContent=`J${s.j0}–J${s.j1}`+(nx&&day>s.mid?` → ${nx.id} (${Math.round(100*(day-s.mid)/(nx.mid-s.mid))} %)`:'');}
// ---- interpolation (projection sur l'homologue du stade suivant)
const mq=[],done=new Set();
function targets(n){const o={};loaded[n.id].traverse(m=>{if(m.isMesh)o[m.name]=m.geometry;});return o;}
function ensureMorph(s,n){if(!MeshBVH||!n||!loaded[n.id]||!root)return;const tg=targets(n);root.traverse(m=>{if(!m.isMesh)return;const id=s.id+'|'+n.id+'|'+m.name;if(done.has(id))return;done.add(id);const t=tg[m.name];if(t)mq.push({m,t,i:0,n:n.id});});}
function morphStep(){return;const j=mq[0];if(!j)return;try{if(!j.t.boundsTree)j.t.boundsTree=new MeshBVH(j.t);const pos=j.m.geometry.attributes.position;if(!j.out)j.out=new Float32Array(pos.count*3);const p=new THREE.Vector3(),tg={point:new THREE.Vector3()};const end=Math.min(pos.count,j.i+12000);
 for(let i=j.i;i<end;i++){p.fromBufferAttribute(pos,i);j.t.boundsTree.closestPointToPoint(p,tg);j.out[3*i]=tg.point.x;j.out[3*i+1]=tg.point.y;j.out[3*i+2]=tg.point.z;}j.i=end;if(end<pos.count)return;
 j.m.geometry.morphAttributes.position=[new THREE.Float32BufferAttribute(j.out,3)];j.m.geometry.morphTargetsRelative=false;j.m.updateMorphTargets();j.m.material.needsUpdate=true;j.m.userData.to=j.n;}catch(e){console.warn(e);}mq.shift();}
function applyMorph(){return; /* interpolation par projection interdite (24/09) */ const s=stageAt(day),n=nextOf(s);if(!root)return;const k=(n&&$('morph').checked)?Math.max(0,Math.min(1,(day-s.mid)/(n.mid-s.mid))):0;root.traverse(m=>{if(m.isMesh&&m.morphTargetInfluences&&m.morphTargetInfluences.length)m.morphTargetInfluences[0]=(m.userData.to===(n&&n.id))?k:0;});if(n&&$('morph').checked&&day>s.mid)ensureMorph(s,n);}
// ---- tubes à topologie commune (session VHE) : maillage fixe (N anneaux × M points), sommets interpolés entre les deux stades voisins
const tubesGroup=new THREE.Group();scene.add(tubesGroup);let TUBES=null;
if(D.tubes){fetch(D.tubes.fichier+'?v='+D.version).then(r=>r.arrayBuffer()).then(buf=>{const all=new Float32Array(buf);window.__tubesData=all;const {N,M,noms,stades}=D.tubes;const idx=[];
 for(let i=0;i<N-1;i++)for(let j=0;j<M;j++){const a=i*M+j,b=i*M+(j+1)%M,c=(i+1)*M+j,d=(i+1)*M+(j+1)%M;idx.push(a,c,b,b,c,d);}
 TUBES=noms.map((n,k)=>{const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.BufferAttribute(new Float32Array(N*M*3),3));g.setIndex(idx);
  const m=new THREE.Mesh(g,new THREE.MeshStandardMaterial({color:new THREE.Color(TUBE_COL(n)),roughness:0.5,metalness:0,side:THREE.DoubleSide}));m.name='tube:'+n;m.visible=false;sysOf[m.name]=TUBE_SYS(n);tubesGroup.add(m);
  return {n,m,base:k*D.tubes.taille,stride:N*M*3,pres:stades.map(s=>!!(D.tubes.presence[n]&&D.tubes.presence[n][s]))};});
 window.__tubes=TUBES;applyVis();}).catch(e=>console.warn('tubes',e));}
function pairAt(d){let a=ST[0],b=null;for(const s of ST){if(s.mid<=d)a=s;else{b=s;break;}}if(!b)return {a,b:a,k:0};if(d<=ST[0].mid)return {a:ST[0],b:ST[0],k:0};return {a,b,k:Math.max(0,Math.min(1,(d-a.mid)/(b.mid-a.mid)))};}
function updateTubes(){if(!TUBES||!D.tubes)return;const on_=$('tubes').checked;const {a,b,k}=pairAt(day);const ia=D.tubes.stades.indexOf(a.id),ib=D.tubes.stades.indexOf(b.id);
 TUBES.forEach(tb=>{const pa=ia>=0&&tb.pres[ia],pb=ib>=0&&tb.pres[ib];let src=null,dst=null,kk=k;
  if(pa&&pb){src=ia;dst=ib;}else if(pa&&k<0.5){src=dst=ia;kk=0;}else if(pb&&k>=0.5){src=dst=ib;kk=0;}
  const vis=on_&&src!==null&&on[sysOf[tb.m.name]]!==false;tb.m.visible=vis;if(!vis)return;
  const key=src+'|'+dst+'|'+kk.toFixed(3);if(tb.key===key)return;tb.key=key;
  const P=tb.m.geometry.attributes.position.array;const o1=tb.base+src*tb.stride,o2=tb.base+dst*tb.stride;const A=window.__tubesData;
  for(let i=0;i<tb.stride;i++)P[i]=A[o1+i]+(A[o2+i]-A[o1+i])*kk;
  tb.m.geometry.attributes.position.needsUpdate=true;tb.m.geometry.computeVertexNormals();tb.m.geometry.computeBoundingSphere();});}
$('tubes').onchange=()=>{applyVis();updateTubes();};if(!D.tubes)$('tubesLbl').style.display='none';
// ---- morphing web : clés à topologie fixe (7 stades) de la session base, sommets interpolés linéairement selon le curseur temps
const morphGroup=new THREE.Group();scene.add(morphGroup);let MORPH=null;
const MORPH_SYS=n=>/^(encephale|moelle)_morph$/.test(n)?'nerveux':/^vertebre_/.test(n)?'colonne':(sysOf[n]||'autres');
const MORPH_COL=(n,c)=>{const h=n==='encephale_morph'?'#4f7fc9':n==='moelle_morph'?'#2f5fb0':/^vertebre_/.test(n)?'#f1ead6':colorOf[n];if(h)return new THREE.Color(h);return Array.isArray(c)?new THREE.Color(c[0],c[1],c[2]):new THREE.Color(c||'#999');};
const morphOn=()=>!!(MORPH&&$('morph').checked);
if(D.morph){$('morphLbl').style.display='';$('autozoom').checked=true;MORPH=D.morph.structures.map(st=>({...st,m:null,keys:null,key:''}));
 MORPH.forEach(st=>{Promise.all([fetch(st.pos+'?v='+D.version).then(r=>r.arrayBuffer()),fetch(st.idx+'?v='+D.version).then(r=>r.arrayBuffer())]).then(([pb,ib])=>{
  st.keys=new Float32Array(pb);const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.BufferAttribute(new Float32Array(st.n*3),3));g.setIndex(new THREE.BufferAttribute(new Uint32Array(ib),1));
  const mat=new THREE.MeshStandardMaterial({color:MORPH_COL(st.nom,st.couleur),roughness:0.55,metalness:0,side:THREE.DoubleSide});
  if(st.alpha<1||D.translucides[st.nom]!=null){mat.transparent=true;mat.opacity=D.translucides[st.nom]!=null?D.translucides[st.nom]:st.alpha;mat.depthWrite=false;}
  st.m=new THREE.Mesh(g,mat);st.m.name='morph:'+st.nom;sysOf[st.m.name]=MORPH_SYS(st.nom);st.m.visible=false;morphGroup.add(st.m);updateMorphWeb(true);}).catch(e=>console.warn('morph',st.nom,e));});}
function updateMorphWeb(force){if(!MORPH)return;const on_=morphOn();const {a,b,k}=pairAt(day);const ia=D.morph.stades.indexOf(a.id),ib=D.morph.stades.indexOf(b.id);
 MORPH.forEach(st=>{if(!st.m)return;const vis=on_&&ia>=0&&ib>=0&&on[sysOf[st.m.name]]!==false;st.m.visible=vis;if(!vis)return;
  const key=ia+'|'+ib+'|'+k.toFixed(3);if(st.key===key&&!force)return;st.key=key;const P=st.m.geometry.attributes.position.array,L=st.n*3,o1=ia*L,o2=ib*L,K=st.keys;
  for(let i=0;i<L;i++)P[i]=K[o1+i]+(K[o2+i]-K[o1+i])*k;st.m.geometry.attributes.position.needsUpdate=true;st.m.geometry.computeVertexNormals();st.m.geometry.computeBoundingSphere();});
 if(root)root.visible=!on_;if(muscGroup)muscGroup.visible=!on_;}
$('morph').onchange=()=>{updateMorphWeb(true);applyVis();};
// ---- temps
$('tl').insertAdjacentHTML('beforeend',ST.map(s=>`<span class="tk" style="left:${(s.mid-26)/(55-26)*100}%">${s.id}</span>`).join(''));
$('day').oninput=()=>{day=+$('day').value;};
$('tab').onclick=()=>document.body.classList.toggle('drawer');canvas.addEventListener('pointerdown',()=>document.body.classList.remove('drawer'));
let fsCss=false;const fsEl=()=>document.fullscreenElement||document.webkitFullscreenElement;
$('fs').onclick=async()=>{if(fsEl()){(document.exitFullscreen||document.webkitExitFullscreen).call(document);return;}
 if(fsCss){fsCss=false;window.parent.postMessage({v3d:'fullscreen',on:false},'*');$('fs').classList.remove('on');return;}
 const el=document.documentElement;const rf=el.requestFullscreen||el.webkitRequestFullscreen;let ok=false;if(rf){try{await rf.call(el);ok=true;}catch(e){}}
 if(!ok&&window.parent!==window){fsCss=true;window.parent.postMessage({v3d:'fullscreen',on:true},'*');$('fs').classList.add('on');}};
document.addEventListener('fullscreenchange',()=>$('fs').classList.toggle('on',!!fsEl()));
$('play').onclick=()=>{playing=!playing;$('play').textContent=playing?'❚❚':'▶';$('play').classList.toggle('on',playing);applyVis();};
function loop(now){requestAnimationFrame(loop);const dt=Math.max(0,Math.min(0.1,(now-last)/1000));last=now;
 const s=stageAt(day);if(!loaded[s.id]&&!morphOn()){load(s,()=>show(s));}else{if(loaded[s.id])show(s);else load(s,()=>show(s));
  if(playing){const nxt=stageAt(day+dt*+$('speed').value);if(loaded[nxt.id]||nxt.id===s.id||morphOn()){day+=dt*+$('speed').value;if(day>ST[ST.length-1].j1+1){if($('loop').checked)day=ST[0].j0-2;else{day=ST[ST.length-1].j1+1;playing=false;$('play').textContent='▶';$('play').classList.remove('on');}}$('day').value=day;}else load(nxt);}
  morphStep();applyMorph();updateTubes();updateMorphWeb();updateInfo();}
 if(fitted&&$('autozoom').checked&&!userMoved&&dt>0){const g=framing(boxAt(day),true);if(!camGoal||camGoal.pos.distanceTo(g.pos)>1e-4||camGoal.tgt.distanceTo(g.tgt)>1e-4)camGoal=g;}
 if(camGoal&&dt>0){sdamp(cam.position,camGoal.pos,camVel.p,0.35,dt);sdamp(ctl.target,camGoal.tgt,camVel.t,0.35,dt);if(cam.position.distanceTo(camGoal.pos)<0.002&&camVel.p.length()<0.002)camGoal=null;}
 $('dayLbl').textContent='J'+day.toFixed(0);ctl.update();R.render(scene,cam);}
window.__dbg=()=>{const b=visibleBox();const sz=b.getSize(new THREE.Vector3());return {pos:cam.position.toArray().map(x=>+x.toFixed(2)),tgt:ctl.target.toArray().map(x=>+x.toFixed(2)),dist:+cam.position.distanceTo(ctl.target).toFixed(2),sz:sz.toArray().map(x=>+x.toFixed(2)),goal:camGoal&&camGoal.pos.toArray().map(x=>+x.toFixed(2)),fov:cam.fov,aspect:+cam.aspect.toFixed(2),cur,userMoved};};
// sondes de test (onglet caché : rAF suspendu) : avancer la boucle à la main
window.__COMMUN=COMMUN;window.__morph=()=>MORPH&&MORPH.map(s=>[s.nom,!!s.m,s.m&&s.m.visible,s.key]);window.__visibles=()=>{const a=[];if(root)root.traverse(o=>{if(o.isMesh&&o.visible)a.push(o.name);});if(muscGroup)muscGroup.traverse(o=>{if(o.isMesh&&o.visible)a.push(o.name);});return a;};window.__ST=ST;window.__setDay=d=>{day=d;$('day').value=d;};window.__step=(n,ms)=>{const raf=window.requestAnimationFrame;window.requestAnimationFrame=()=>0;let t=last;for(let i=0;i<n;i++){t+=ms;loop(t);}window.requestAnimationFrame=raf;};
load(ST[0],()=>{show(ST[0]);requestAnimationFrame(loop);});
</script></body></html>
"""
out_dir = os.path.join("embryons_3D", "site_3dht" if EMBED else "site_3dh")
os.makedirs(out_dir, exist_ok=True)


def ecrire_tubes(out_dir):
    """Tubes à topologie commune (session VHE) : embryons_3D/tubes_morph.npz {structure: (7 stades, N anneaux, M points, xyz mm)}
    → out_dir/tubes.bin (float32, structures concaténées) + métadonnées dans la page. None si absent."""
    npz, js = os.path.join("embryons_3D", "tubes_morph.npz"), os.path.join("embryons_3D", "tubes_morph.json")
    if not (os.path.exists(npz) and os.path.exists(js)):
        return None
    import numpy as np
    z = np.load(npz); j = json.load(open(js, encoding="utf-8"))
    noms = list(z.files)
    n_st, N, M = z[noms[0]].shape[:3]
    arr = np.concatenate([z[n].astype("<f4").reshape(1, -1) for n in noms], 0)
    with open(os.path.join(out_dir, "tubes.bin"), "wb") as f:
        f.write(arr.tobytes())
    return {"noms": noms, "stades": j["stades"], "N": int(N), "M": int(M), "presence": j["presence"], "fichier": "tubes.bin",
            "taille": int(n_st * N * M * 3)}


data["tubes"] = ecrire_tubes(out_dir)


def ecrire_morph(out_dir):
    """Clés de morphing de la session base : embryons_3D/morph_web/morph_web.json + <structure>.npz {positions (7,N,3) mm, faces (F,3)}
    → out_dir/morph/<structure>_pos.bin (float32) + <structure>_idx.bin (uint32) ; métadonnées dans la page. None si absent."""
    src = os.environ.get("MORPH_DIR", os.path.join("embryons_3D", "morph_web"))
    js = os.path.join(src, "morph_web.json")
    if not os.path.exists(js):
        return None
    import numpy as np
    j = json.load(open(js, encoding="utf-8"))
    dst = os.path.join(out_dir, "morph"); os.makedirs(dst, exist_ok=True)
    structs = []
    for st in j.get("structures", []):
        f = os.path.join(src, st.get("fichier") or (st["nom"] + ".npz"))
        if not os.path.exists(f):
            continue
        z = np.load(f); pos = np.ascontiguousarray(z["positions"], dtype="<f4"); fac = np.ascontiguousarray(z["faces"], dtype="<u4")
        n_st, n, _ = pos.shape
        with open(os.path.join(dst, st["nom"] + "_pos.bin"), "wb") as fh:
            fh.write(pos.tobytes())
        with open(os.path.join(dst, st["nom"] + "_idx.bin"), "wb") as fh:
            fh.write(fac.tobytes())
        structs.append({"nom": st["nom"], "pos": f"morph/{st['nom']}_pos.bin", "idx": f"morph/{st['nom']}_idx.bin", "n": int(n), "f": int(len(fac)),
                        "n_st": int(n_st), "couleur": st.get("couleur"), "alpha": st.get("alpha", 1.0)})
    if not structs:
        return None
    return {"stades": j["stades"], "structures": structs}


data["morph"] = ecrire_morph(out_dir)
fn = "viewer.html" if EMBED else ("index_local.html" if LOCAL else "index.html")
with open(os.path.join(out_dir, fn), "w", encoding="utf-8", newline="") as f:
    f.write(HTML.replace("__DATA__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/")))
if not LOCAL:
    with open(os.path.join(out_dir, VERSION_FILE), "w") as f:
        f.write(data["version"])
print("OK", os.path.join(out_dir, fn), "version", data["version"], "| stades", [s["id"] for s in stades])
