# -*- coding: utf-8 -*-
"""Publication automatique de harcelon.fr/3dht au fil de l'avancée.

    python embryo3d/publier_auto.py            # surveille en continu (toutes les 30 s)
    python embryo3d/publier_auto.py --une-fois # une passe puis sort

Surveille embryons_3D/modeles/ (manifests et maillages), embryo3d/agregateur.py et embryons_3D/taches_etat.json.
À chaque changement stable (plus rien ne bouge depuis 20 s, pour ne pas publier une livraison à moitié écrite) :
agregateur.py --site → envoi des dossiers modeles/ modifiés → index.html et version.txt EN DERNIER (les pages ouvertes
se rechargent alors d'elles-mêmes). Un modèle « publiable": false » n'est jamais envoyé (filtré par agregateur.py) ;
les dossiers retirés du paquet sont aussi retirés du serveur.
"""
import hashlib
import os
import subprocess
import sys
import time
from datetime import datetime

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)
SITE = os.path.join("embryons_3D", "site_3dht")
DISTANT = "harcelon:/home/ubuntu/havre-vps/site/3dht"
HOTE, CHEMIN = DISTANT.split(":", 1)
SURVEILLES = [os.path.join("embryons_3D", "modeles"), os.path.join("embryo3d", "agregateur.py"), os.path.join("embryons_3D", "taches_etat.json")]
PERIODE, CALME = 30, 20
LOG = os.path.join("embryons_3D", "publier_auto.log")
ETAT_PUB = os.path.join("embryons_3D", "publier_auto_etat.json")   # signatures des dossiers réellement envoyés au serveur


def lire_envoyes():
    import json
    try:
        return json.load(open(ETAT_PUB, encoding="utf-8"))
    except Exception:
        return {}


def ecrire_envoyes(d):
    import json
    json.dump(d, open(ETAT_PUB, "w", encoding="utf-8"), indent=1)


def journal(t):
    ligne = f"{datetime.now():%d/%m %H:%M:%S} {t}"
    print(ligne, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(ligne + "\n")


def empreinte():
    """(signature, date du dernier changement) des fichiers surveillés."""
    items, dernier = [], 0.0
    for chemin in SURVEILLES:
        fichiers = []
        if os.path.isdir(chemin):
            for r, _, fs in os.walk(chemin):
                fichiers += [os.path.join(r, f) for f in fs if not f.startswith("verif_")]
        elif os.path.exists(chemin):
            fichiers = [chemin]
        for f in fichiers:
            st = os.stat(f)
            items.append(f"{f}|{st.st_size}|{int(st.st_mtime)}")
            dernier = max(dernier, st.st_mtime)
    return hashlib.md5("\n".join(sorted(items)).encode()).hexdigest(), dernier


def signatures_site():
    """signature par dossier de modèle dans le paquet site."""
    base = os.path.join(SITE, "modeles"); out = {}
    if os.path.isdir(base):
        for n in os.listdir(base):
            d = os.path.join(base, n)
            out[n] = hashlib.md5("\n".join(sorted(f"{f}|{os.path.getsize(os.path.join(d, f))}|{int(os.path.getmtime(os.path.join(d, f)))}"
                                                   for f in os.listdir(d))).encode()).hexdigest()
    return out


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding="utf-8", errors="replace",
                       env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    if r.returncode:
        raise RuntimeError(f"{cmd} → {r.returncode} {r.stderr.strip()[-300:]}")
    return r.stdout


def distants():
    """{chemin relatif sous modeles/ : (taille, date)} des fichiers en ligne."""
    out = run("ssh " + HOTE + " \"cd " + CHEMIN + " && mkdir -p modeles && cd modeles && find . -type f -printf '%P %s %T@\\n'\"")
    d = {}
    for ligne in out.splitlines():
        if ligne.count(" ") >= 2:
            chemin, taille, date = ligne.rsplit(" ", 2)
            d[chemin] = (int(taille), float(date))
    return d


def publier(avant=None):
    """Synchronisation fichier par fichier (taille comparée au serveur) : on n'efface plus un dossier avant de le renvoyer,
    un stade en ligne n'est donc jamais incomplet ; les nouveaux fichiers passent avant la page, les fichiers retirés après."""
    run(f'python embryo3d/agregateur.py --site "{SITE}"')
    base = os.path.join(SITE, "modeles")
    locaux = {}
    for r, _, fs in os.walk(base):
        for f in fs:
            rel = os.path.relpath(os.path.join(r, f), base).replace("\\", "/")
            st = os.stat(os.path.join(r, f))
            locaux[rel] = (st.st_size, st.st_mtime)
    enligne = distants()
    # à envoyer : absent en ligne, taille différente, ou régénéré après l'envoi (même taille, contenu neuf : PLY à nombre de sommets fixe)
    a_envoyer = sorted(k for k, (t, m) in locaux.items() if k not in enligne or enligne[k][0] != t or m > enligne[k][1] + 2)
    par_dossier = {}
    for k in a_envoyer:
        par_dossier.setdefault(k.split("/")[0], []).append(k)
    for dossier, fichiers in par_dossier.items():
        run(f'ssh {HOTE} "mkdir -p {CHEMIN}/modeles/{dossier}"')
        for i in range(0, len(fichiers), 40):                     # scp groupé par lots de 40 fichiers
            lot = " ".join(f'"{base}/{k}"' for k in fichiers[i:i + 40])
            run(f'scp -q {lot} {DISTANT}/modeles/{dossier}/')
    if os.path.exists(os.path.join("embryons_3D", "taches_etat.json")):
        run(f'scp -q embryons_3D/taches_etat.json {DISTANT}/taches_etat.json')
    run(f'scp -q "{SITE}/index.html" "{SITE}/version.txt" {DISTANT}/')      # la page, une fois tous ses fichiers en place
    en_trop = sorted(k for k in enligne if k not in locaux)
    for i in range(0, len(en_trop), 60):
        run(f'ssh {HOTE} "cd {CHEMIN}/modeles && rm -f ' + " ".join(f"'{k}'" for k in en_trop[i:i + 60]) + '"')
    run(f'ssh {HOTE} "cd {CHEMIN}/modeles && find . -type d -empty -delete"')
    v = open(os.path.join(SITE, "version.txt")).read().strip()
    journal(f"publié {v}" + (f" · {len(a_envoyer)} fichiers envoyés ({', '.join(par_dossier)})" if a_envoyer else "")
            + (f" · {len(en_trop)} retirés" if en_trop else ""))
    return {}


def main():
    une_fois = "--une-fois" in sys.argv
    publie, _ = (None, 0) if une_fois else empreinte()
    site = {}
    journal("surveillance démarrée" if not une_fois else "passe unique")
    while True:
        sig, dernier = empreinte()
        if sig != publie and (une_fois or time.time() - dernier >= CALME):
            try:
                site = publier(site)
                publie = sig
            except Exception as e:
                journal(f"ÉCHEC : {e}")
        if une_fois:
            break
        time.sleep(PERIODE)


if __name__ == "__main__":
    main()
