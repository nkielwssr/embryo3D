# -*- coding: utf-8 -*-
"""
Lanceur local des tâches de l'agrégateur (embryons_3D/agregateur.html).

    python embryo3d/lanceur.py [port=8791]

La page (servie en localhost) appelle ce serveur pour :
  GET  /etat                → état de toutes les tâches (embryons_3D/taches_etat.json)
  POST /lancer  {id}        → exécute la commande de la tâche (agregateur.json) dans un sous-processus,
                              journal dans embryons_3D/taches_logs/<id>.log ; statut en_cours → termine / echec
  POST /confier {id, instance}  → attribue la tâche à une instance (session Claude nommée, personne…)
  POST /statut  {id, statut}    → a_faire | en_cours | termine | echec  (marquage manuel)
  GET  /log/<id>            → journal texte de la tâche
Fichier d'état partagé : les autres sessions Claude peuvent l'éditer directement (même format).
Option --publier [cible scp] : à chaque changement du fichier d'état (y compris par une autre session), copie
vers le site (défaut harcelon:/home/ubuntu/havre-vps/site/3dht/taches_etat.json) avec l'horodatage "_publie" ;
la page en ligne le relit toutes les 30 s.
"""
import json
import os
import subprocess
import sys
import threading
import time
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)
ETAT = os.path.join("embryons_3D", "taches_etat.json")
TACHES = os.path.join("embryons_3D", "agregateur.json")
LOGS = os.path.join("embryons_3D", "taches_logs")
os.makedirs(LOGS, exist_ok=True)
verrou = threading.Lock()
procs = {}


def lire_etat():
    if os.path.exists(ETAT):
        try:
            return json.load(open(ETAT, encoding="utf-8"))
        except Exception:
            pass
    return {}


def ecrire_etat(e):
    with open(ETAT, "w", encoding="utf-8") as f:
        json.dump(e, f, ensure_ascii=False, indent=1)


def maj(tid, **champs):
    with verrou:
        e = lire_etat()
        t = e.get(tid, {})
        t.update(champs)
        t["maj"] = datetime.now().strftime("%d/%m/%Y %H:%M")
        e[tid] = t
        ecrire_etat(e)
        return t


def taches():
    if not os.path.exists(TACHES):
        return {}
    d = json.load(open(TACHES, encoding="utf-8"))
    return {t["id"]: t for t in d.get("todo", []) if "id" in t}


def lancer(tid, cmd):
    log = os.path.join(LOGS, f"{tid}.log")
    f = open(log, "a", encoding="utf-8")
    f.write(f"\n=== {datetime.now():%d/%m/%Y %H:%M:%S} $ {cmd}\n"); f.flush()
    p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, cwd=RACINE,
                         env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    procs[tid] = p
    maj(tid, statut="en_cours", instance=f"lanceur local (PID {p.pid})", cmd=cmd, lance=datetime.now().strftime("%d/%m/%Y %H:%M"), log=log)

    def attendre():
        code = p.wait()
        f.write(f"=== fin, code {code}\n"); f.close()
        maj(tid, statut="termine" if code == 0 else "echec", code=code)
        procs.pop(tid, None)
    threading.Thread(target=attendre, daemon=True).start()


class H(BaseHTTPRequestHandler):
    def _rep(self, code, corps, ctype="application/json"):
        data = corps if isinstance(corps, bytes) else json.dumps(corps, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self._rep(204, b"")

    def do_GET(self):
        if self.path.startswith("/etat"):
            e = lire_etat()
            for tid, p in list(procs.items()):
                if tid in e:
                    e[tid]["vivant"] = p.poll() is None
            return self._rep(200, {"etat": e, "lanceur": True, "heure": datetime.now().strftime("%H:%M:%S")})
        if self.path.startswith("/log/"):
            tid = self.path[5:].split("?")[0]
            log = os.path.join(LOGS, f"{tid}.log")
            if os.path.exists(log):
                return self._rep(200, open(log, "rb").read()[-20000:], "text/plain")
            return self._rep(404, b"pas de journal", "text/plain")
        return self._rep(404, {"erreur": "inconnu"})

    def do_POST(self):
        n = int(self.headers.get("Content-Length", 0))
        try:
            corps = json.loads(self.rfile.read(n).decode("utf-8") or "{}")
        except Exception:
            corps = {}
        tid = str(corps.get("id", ""))
        if not tid:
            return self._rep(400, {"erreur": "id manquant"})
        if self.path.startswith("/lancer"):
            t = taches().get(tid)
            cmd = corps.get("cmd") or (t or {}).get("cmd")
            if not cmd:
                return self._rep(400, {"erreur": "pas de commande pour cette tâche : la confier à une instance"})
            if tid in procs and procs[tid].poll() is None:
                return self._rep(409, {"erreur": "déjà en cours"})
            lancer(tid, cmd)
            return self._rep(200, {"ok": True, "etat": lire_etat().get(tid)})
        if self.path.startswith("/confier"):
            inst = str(corps.get("instance", "")).strip()
            return self._rep(200, {"ok": True, "etat": maj(tid, instance=inst, statut=corps.get("statut", "en_cours"))})
        if self.path.startswith("/statut"):
            st = corps.get("statut", "a_faire")
            if st not in ("a_faire", "en_cours", "termine", "echec"):
                return self._rep(400, {"erreur": "statut inconnu"})
            return self._rep(200, {"ok": True, "etat": maj(tid, statut=st)})
        return self._rep(404, {"erreur": "inconnu"})

    def log_message(self, *a):
        pass


CIBLE_PUBLICATION = "harcelon:/home/ubuntu/havre-vps/site/3dht/taches_etat.json"


def publier_en_boucle(cible):
    """Surveille le fichier d'état et le pousse (scp) dès qu'il change."""
    dernier = None
    tmp = os.path.join("embryons_3D", ".taches_etat_publie.json")
    while True:
        try:
            m = os.path.getmtime(ETAT) if os.path.exists(ETAT) else None
            if m and m != dernier:
                e = lire_etat()
                e["_publie"] = datetime.now().strftime("%d/%m/%Y %H:%M")
                with open(tmp, "w", encoding="utf-8") as f:
                    json.dump(e, f, ensure_ascii=False, indent=1)
                r = subprocess.run(["scp", "-q", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", tmp, cible],
                                   capture_output=True, text=True, timeout=120)
                print(f"[{datetime.now():%H:%M:%S}] publication de l'état -> {cible} : "
                      f"{'ok' if r.returncode == 0 else 'echec ' + r.stderr.strip()}", flush=True)
                if r.returncode == 0:
                    dernier = m
        except Exception as ex:
            print("publication :", ex, flush=True)
        time.sleep(10)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--publier"]
    port = int(args[0]) if args and args[0].isdigit() else 8791
    if "--publier" in sys.argv:
        cible = args[1] if len(args) > 1 else CIBLE_PUBLICATION
        threading.Thread(target=publier_en_boucle, args=(cible,), daemon=True).start()
    print(f"lanceur des tâches : http://127.0.0.1:{port}/etat  (état : {ETAT})")
    ThreadingHTTPServer(("127.0.0.1", port), H).serve_forever()
