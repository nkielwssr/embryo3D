# -*- coding: utf-8 -*-
"""Test de la table hikspoors_nomenclature.json sur les noms réels des PDF Hikspoors (inventaire du relais local, 26/09 : CS9, CS13,
CS18_NEWvalves, CS20_NEWvalves, CS23_NEWvalves, plus les variantes des versions d'origine). Sans données, quelques secondes.

    python embryo3d/test_hikspoors_nomenclature.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hikspoors_modele as hm                                   # noqa: E402

IGNORER = None
ATTENDU = {
    # CS9
    "arterial_plexus": "plexus_arteriel", "cardiac_jelly": "gelee_cardiaque", "coelom": "coelome",
    "gut_yolk_sac": "intestin_et_vesicule_vitelline", "myocard_LV": "myocarde_ventricule_gauche", "neural_plate": "tube_neural",
    "pericard": "pericarde", "scale_cube_200um": IGNORER, "somites_somitomeres": "somites", "transverse_septum": "septum_transversum",
    "venous_plexus": "plexus_veineux",
    # CS13
    "DMP": "protrusion_mesenchymateuse_dorsale", "L_HCC": "canal_hepatocardiaque_gauche", "R_HCC": "canal_hepatocardiaque_droit",
    "L_Umb_vein": "veine_ombilicale_gauche", "R_Umb_vein": "veine_ombilicale_droit", "L_Vit_vein": "veine_vitelline_gauche",
    "R_Vit_vein": "veine_vitelline_droit", "L_cardinal_vein": "veine_cardinale_gauche", "R_cardinal_vein": "veine_cardinale_droit",
    "L_common_card_vein": "veine_cardinale_commune_gauche", "R_common_card_vein": "veine_cardinale_commune_droit",
    "OFT_cardiac_jelly": "gelee_cardiaque_voie_efferente", "PAAs_2": "arc_aortique_2", "PAAs_3": "arc_aortique_3", "PAAs_4": "arc_aortique_4",
    "SAN": "noeud_sinusal", "aortic_sac": "sac_aortique", "dorsal_aorta": "aorte_dorsale", "epicard": "epicarde",
    "inf_endocardial_cushion": "coussin_endocardique_inferieur", "sup_endocardial_cushion": "coussin_endocardique_superieur",
    "liver": "foie", "loop_wire_heart_tube": "axe_tube_cardiaque", "lumen_LA": "cavite_oreillette_gauche", "lumen_LV": "cavite_ventricule_gauche",
    "lumen_OFT": "cavite_voie_efferente", "lumen_RA": "cavite_oreillette_droite", "lumen_RV": "cavite_ventricule_droit",
    "lumen_venous_sinus": "cavite_sinus_veineux", "myocard_AV_canal": "myocarde_canal_atrioventriculaire", "myocard_OFT": "myocarde_voie_efferente",
    "myocard_RV": "myocarde_ventricule_droit", "myocard_atriums": "myocarde_oreillettes", "myocard_venous_sinus": "myocarde_sinus_veineux",
    "neural_tube": "tube_neural", "pericardial_reflection": "reflexion_pericardique", "pulmonary_vein": "veine_pulmonaire", "somites": "somites",
    # CS18 / CS20 / CS23 (NEWvalves et originales)
    "Ao_non_adj_leaflet": "valve_aortique_non_adj", "Ao_parietal_leaflet": "valve_aortique_parietal", "Ao_septal_leaflet": "valve_aortique_septal",
    "pulm_non_adj_leaflet": "valve_pulmonaire_non_adj", "pulm_parietal_leaflet": "valve_pulmonaire_parietal",
    "pulm_septal_leaflet": "valve_pulmonaire_septal", "CCS": "systeme_conduction", "CCS_1": "systeme_conduction", "CSS": "systeme_conduction",
    "L_spinal_ganglia": "ganglions_spinaux_gauche", "R_spinal_ganglia": "ganglions_spinaux_droit", "NCCs": "cellules_crete_neurale",
    "PAAs_6": "arc_aortique_6", "asc_Ao_wall": "paroi_aorte_ascendante", "Asc_Ao_wall": "paroi_aorte_ascendante",
    "ascending_aorta": "aorte_ascendante", "cardinal_veins": "veine_cardinale", "common_cardinal_veins": "veine_cardinale_commune",
    "left_coronary_artery": "arteres_coronaires", "coronary_arteries": "arteres_coronaires",
    "lumen_OFT_subaortic_part": "cavite_voie_efferente_sous_aortique", "lumen_OFT_subpulm_part": "cavite_voie_efferente_sous_pulmonaire",
    "lumen_pulmonary_vein": "veine_pulmonaire", "myocard_pulmonary_veins": "myocarde_veines_pulmonaires", "lungs": "poumons",
    "main_branches_lungs": "poumons", "lung_lobes": "lobes_pulmonaires", "mitral_valve": "valve_mitrale", "tricuspid_valve": "valve_tricuspide",
    "musc_ventricular_septum": "septum_interventriculaire", "myocard_outlet_septum": "septum_voie_efferente",
    "myocardial_DMP": "myocarde_protrusion_mesenchymateuse_dorsale", "primary_atrial_septum": "septum_primum",
    "secondary_atrial_septum": "septum_secundum", "pulm_trunk_wall": "paroi_tronc_pulmonaire", "pulmonary_arteries": "arteres_pulmonaires",
    "pulmonary_trunk": "tronc_pulmonaire", "venous_valves": "valvules_sinusales", "gut": "intestin", "azygos_venous_system": "veine_azygos",
    "cranial_veins": "veines_craniennes", "inf_caval_vein": "veine_cave_inferieure", "sup_caval_vein": "veine_cave_superieure",
    "lumen_coronary_sinus": "sinus_coronaire",
    # noms du test synthétique (forme « à la Hikspoors » supposée avant l'inventaire) : ne pas régresser
    "Myocardium_LV": "myocarde_ventricule_gauche", "Lumen_RV": "cavite_ventricule_droit", "OFT_myocardium": "myocarde_voie_efferente",
    "Sinus_venosus": "coeur_sinus_veineux", "Umbilical_vein_L": "veine_ombilicale_gauche", "Dorsal_aorta_right": "aorte_dorsale_droit",
    "Aortic_arch_3_left": "arc_aortique_3_gauche", "Septum_transversum": "septum_transversum", "Liver": "foie",
    # pièges des anciens motifs (\w* traverse '_', « ear » dans heart, « card » dans cardinal)
    "heart_tube": "coeur_tube_cardiaque", "myocardium": "myocarde_coeur", "atrial_septum": "septum_interatrial",
    # retour 2 du relais (26/09 soir) : CS10-CS12 et CS14-CS17
    "lumen_IFT": "cavite_voie_afferente", "myocard_IFT": "myocarde_voie_afferente", "IFT": "coeur_voie_afferente",
    "Ao_swelling": "bourrelet_valvulaire_aortique", "pulm_swelling": "bourrelet_valvulaire_pulmonaire", "scale_cube___200um": IGNORER,
}
SYSTEMES = {"pericard": "annexes", "epicard": "vasculaire", "L_spinal_ganglia": "nerveux", "NCCs": "tissus", "coelom": "annexes",
            "lung_lobes": "respiratoire", "R_cardinal_vein": "vasculaire"}
GROUPES = {"R_cardinal_vein": "veines", "L_common_card_vein": "veines", "pericard": None, "epicard": "coeur", "PAAs_3": "arteres",
           "lumen_venous_sinus": "cavites_cardiaques", "myocard_venous_sinus": "coeur", "loop_wire_heart_tube": None,
           "OFT_cardiac_jelly": "coeur", "L_Vit_vein": "veines", "lumen_OFT_subpulm_part": "cavites_cardiaques", "gut": "intestin"}


def main():
    nomen = hm.Nomenclature()
    erreurs = []
    for src, att in ATTENDU.items():
        r = nomen.chercher(src)
        obtenu = None if r.get("ignorer") else r["nom"]
        if obtenu != att:
            erreurs.append("%-28s -> %-40s attendu %s   (règle %s)" % (src, obtenu, att, r.get("regle")))
        elif att is not None and r.get("non_reconnu"):
            erreurs.append("%s : reconnu par défaut seulement" % src)
    for src, sy in SYSTEMES.items():
        if nomen.chercher(src)["systeme"] != sy:
            erreurs.append("%s : système %s, attendu %s" % (src, nomen.chercher(src)["systeme"], sy))
    for src, g in GROUPES.items():
        if nomen.chercher(src)["groupe"] != g:
            erreurs.append("%s : groupe %s, attendu %s" % (src, nomen.chercher(src)["groupe"], g))
    # par couleur (surfaces texturées scindées)
    for src, rgb, att in (("gut_yolk_sac", [174, 174, 174], "intestin"), ("gut_yolk_sac", [120, 119, 119], "vesicule_vitelline"),
                          ("myocard_LV", [127, 173, 87], "courbure_interne_lv")):
        r = nomen.par_couleur_regle(src, rgb)
        if r is None or r["nom"] != att:
            erreurs.append("par couleur %s %s -> %s, attendu %s" % (src, rgb, r and r["nom"], att))
    print("%d noms, %d erreurs" % (len(ATTENDU), len(erreurs)))
    if erreurs:
        print("ÉCHEC :\n  " + "\n  ".join(erreurs)); sys.exit(1)
    print("OK : nomenclature Hikspoors conforme à l'inventaire du relais")


if __name__ == "__main__":
    main()
