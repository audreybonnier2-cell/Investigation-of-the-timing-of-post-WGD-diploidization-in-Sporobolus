#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
from ete3 import Tree
# SYS => Gestion des arguments de la ligne de commande
# OS => Gestion des chemins de fichiers et dossiers
# ETE3 Tree => Bibliothèque pour lire et manipuler les arbres Newick
# Espèces considérées Smaritimus et Salterniflorus, script valable pour 2sp
# Notes perso : SP et BM2 vérifier conjointement

def main():
    # Vérification des arguments (Entrée/Sortie)
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_dir> <output_csv>")
        return # Arrêt script si arguments manquent
    input_dir = sys.argv[1] # Dossier des fichiers .nwk
    output_csv = sys.argv[2] # Nom du fichier CSV résultat
    # Ouverture fichier CSV ("w" pour write)
    with open(output_csv, "w") as f_out:

# Écriture entête CSV
        header = [
            "OG", "Bloc", "Smaritimus", "Salterniflorus", "Sum",
            "Total_Dup", "Total_Spec", "Total_T", "Root_Type", "Dist_1st_Spec", "Dist_2nd_Spec", 
            "nb_D_pre1st_S", "nb_D_post1st_S", "1SD_sp1", "1SD_sp2",
            "nb_D_pre2nd_S", "nb_D_post2nd_S", "2SD_sp1", "2SD_sp2",
            "nb_D_pre3rd_S", "nb_D_post3rd_S", "3SD_sp1", "3SD_sp2",
            "nb_D_pre4th_S", "nb_D_post4th_S", "4SD_sp1", "4SD_sp2",
            "S1", "S2", "S3", "S4", 
            "type_scenario", 
            "total_bloc_DFND", "total_bloc_DFD", "total_bloc_absD", "total_bloc_autres", 
            "Newick"
        ]
        # séparation CSV
        f_out.write(",".join(header) + "\n")
        # Initialisation des 3 compteurs pour analyse complète
        total_trees_scanned = 0
        total_ogs_with_blocs = 0
        total_blocs_count = 0
        # Boucle de parcours all fichiers du dossier d'entrée, - ordre alphabétique
        for filename in sorted(os.listdir(input_dir)):
            # Filtre sécurité ignore fichiers - non arbres au format texte/newick
            if not filename.endswith((".newick", ".nwk", ".tree", ".txt")):
                # Si fichier mauvaise exten, passe au fichier suivant
                continue
                
            # Décompte du compteur d'arbres fichier valide
            total_trees_scanned += 1
            # Reconstr chemin d'accès absolu vers le fichier d'arbre actuel
            nwk_path = os.path.join(input_dir, filename)
            # Extraction id OG + retrait extension du fichier
            og_id = os.path.splitext(filename)[0]

# Bloc de sécurité (try) 
#empêche crash si fichier Newick malformé
            try:
                # Chargement arbre (Format 1 pour lire les labels S et D)
                t = Tree(nwk_path, format=1)
                # Compteur blocs pour cet OG
                bloc_id = 1
                # og_written => pour savoir si l'OG a écrit au moins 1 bloc
                og_written = False
##################################################################################################
# ÉTAPE 1 : PARCOURS TOP-DOWN (de la racine vers les feuilles)
##################################################################################################

                for node in t.traverse(strategy="preorder"):
                    # SÉCURITÉ : Si parent déjà validé => ignore ses enfants (évite doublons)
                    if hasattr(node, "stop"):
                        # On saute le nœud et passe au suivant dans le parcours
                        continue
                    # Sécurité : On ignore les feuilles de l'arbre car un bloc est défini par un nœud interne (un ancêtre)
                    if node.is_leaf():
                        # Passe au nœud suivant
                        continue
##################################################################################################
# ÉTAPE 2 : TEST DE MONOPHYLIE LOCALE
##################################################################################################
# Création ensemble (set) vide pour lister les types d'espèces sous ce nœud sans doublons
                    species_found = set()
                    # Boucle pour analyser individuellement chaque feuille (gène) descendant du nœud actuel
                    for l in node.get_leaves():
                        # Conversion nom feuille (évite problèmes de casse)
                        lname = l.name.lower()
                        # Si mot "smaritimus" détecté @nom de la feuille
                        if "smaritimus" in lname:
                            # alors ajout tag "smar" @ensemble d'espèces
                            species_found.add("smar")
                        # Si "salterniflorus" détecté @nom de la feuille
                        elif "salterniflorus" in lname:
                            # alors ajout tag "salt" @ensemble d'espèces
                            species_found.add("salt")
                        # Si autre plante (Acomosus, Zmays, etc) => intrus phylogénétique
                        else:
                            species_found.add("intrus") # Marque la présence d'une autre espèce
# CONDITION : Filtre strictement les deux espèces ET aucun intrus (Le plus gros bloc pur)
                    if "smar" in species_found and "salt" in species_found and "intrus" not in species_found:
                        
                        #  NOUVEAU : TEST DE TRANSFERT 
                        # vérif si  nœud "T" existe dans ce bloc
                        has_transfer = any(n.name and n.name.upper().startswith("T") for n in node.traverse())
                        # Si un transfert est trouvé, on ignore ce bloc (continue)
                        if has_transfer:
                            continue
                                    
                        #ajout down19
                        all_S = [n for n in node.traverse() if not n.is_leaf() and n.name and n.name.upper().startswith("S")]
                        total_spec_bloc = len(all_S)
                        
# Ajout down14+19 => blocs avec 0 spéciation sont ignorés (selon scénario établi)
                        if total_spec_bloc == 0:
                            continue
                            
                        #variable sécurité pour root-type down19
                        root_bloc = node

# VERROU : (Pour n'avoir que Bloc 1 et Bloc 4 ; évite les sous blocs redondants)
                        # On marque ce noeud ET TOUS ses descendants pour qu'ils ne soient plus analysés
                        for descendant in node.traverse():
                            # marquage tous les nœuds enfant avec attribut "stop" =>bloquer @parcours principal
                            descendant.add_feature("stop", True)
                        # len([...]) => Compte nb de gènes individuels de chaque espèce dans ce bloc isolé
                        # Comptage nb sequences Smaritimus @bloc
                        smar_count = len([l for l in node.get_leaves() if "smaritimus" in l.name.lower()])
                        # Comptage nb sequences Salterniflorus @bloc
                        salter_count = len([l for l in node.get_leaves() if "salterniflorus" in l.name.lower()])
                        # Calcul somme totale sequences cibles réunis @bloc (colonne "Sum") à verfier changement "gene_sum"
                        seq_sum = smar_count + salter_count
##################################################################################################
# ETAPE 3 COMPTAGES GLOBAUX DU BLOC
##################################################################################################
# Extraction et comptage duplications (noeuds internes dont le nom commence par 'D')
                        total_dup_bloc = sum(1 for n in node.traverse() if not n.is_leaf() and n.name and n.name.upper().startswith("D"))
                        # COMPTEUR T : extraction et décompte des nœuds internes commençant par "T" => ajout down18
                        total_t_bloc = sum(1 for n in node.traverse() if not n.is_leaf() and n.name and n.name.upper().startswith("T"))
                        
##################################################################################################
# ETAPE 4 Identification Root_Type
##################################################################################################
# Vérification si le nœud racine possède un nom de notation
                        if root_bloc.name:
                            # Si nom commence par D (Duplication)
                            if root_bloc.name.upper().startswith("D"):
                                # racine est officiellement classée comme "D"
                                root_type = "D"
                            # Si nom commence par S (Spéciation)
                            elif root_bloc.name.upper().startswith("S"):
                                # racine est officiellement classée comme "S"
                                root_type = "S"
                            elif root_bloc.name.upper().startswith("T"):
                                root_type = "T"
                            # Si notation est autre ou absente
                            else:
                                # On recopie la valeur textuelle trouvée par sécurité
                                root_type = root_bloc.name
                        else:
                            # Si nœud anonyme, =>NA
                            root_type = "N/A"
##################################################################################################
# ÉTAPE 5-6 : CLASSIFICATION C - TRI PAR DISTANCE ET CALCULS DE DUPLICATIONS CIBLÉES (PRE / POST SPÉCIATION)
##################################################################################################
# DEPLACEMENT ICI DOWN22: Fonction locale pour mesurer la profondeur topologique (en nombre de nœuds intermédiaires)
                        def get_topological_depth(target, root_node):
                            depth = 0
                            curr = target
                            while curr and curr != root_node:
                                depth += 1
                                curr = curr.up
                            return depth
                            
                        # Tri de toutes les spéciations vues par distance évolutive croissante mis à jour down18 => profondeur topo
                        all_S_sorted = sorted(all_S, key=lambda x: get_topological_depth(x, root_bloc))
                        
                        # modif down21 Extraction des distances pour S1 (garantie par le filtre) et S2 (conditionnelle)
                        if len(all_S_sorted) >= 2:
                            dist_1st_spec = all_S_sorted[0].get_distance(root_bloc) # Index 0 => 1st speciation
                            dist_2nd_spec = all_S_sorted[1].get_distance(root_bloc) # Index 1 => 2nd speciation
                        else:
                            dist_1st_spec = all_S_sorted[0].get_distance(root_bloc) if len(all_S_sorted) >= 1 else 0.0
                            dist_2nd_spec = 0.0    
                        # Listes stockage temporaires 
                        pre_counts = []
                        post_counts = []
                        s_classifications = []
                        SD_sp1_list = []
                        SD_sp2_list = []  
                        
                        # ~~ down23 : get_robust_dup_distances
                        # Logique :
                        # 1) Vérifie que sp1 ET sp2 sont présents sous spec_node → sinon N/A, N/A
                        # 2) Cherche la 1ère dup dont les feuilles contiennent UNIQUEMENT sp1 (côté sp1 pur)
                        # 3) Cherche la 1ère dup dont les feuilles contiennent UNIQUEMENT sp2 (côté sp2 pur)
                        # 4) Si les DEUX existent → retourne (dist_sp1, dist_sp2) via node.dist
                        #    Si l'une des deux manque → retourne ("N/A", "N/A")
                        def get_robust_dup_distances(spec_node, sp1_name, sp2_name):
                            leaves_spec = spec_node.get_leaves()
                            has_sp1 = any(sp1_name in l.name.lower() for l in leaves_spec)
                            has_sp2 = any(sp2_name in l.name.lower() for l in leaves_spec)
                            if not (has_sp1 and has_sp2):
                                return "N/A", "N/A"
                        
                            def first_specific_dup(start_node, target_sp, other_sp):
                                # BFS depuis les enfants de start_node
                                # Cherche le 1er nœud D dont les feuilles ont target_sp mais PAS other_sp
                                queue_d = list(start_node.children)
                                while queue_d:
                                    c = queue_d.pop(0)
                                    if c.is_leaf():
                                        continue
                                    if c.name and c.name.upper().startswith("D"):
                                        c_leaves = c.get_leaves()
                                        only_target = (
                                            any(target_sp in l.name.lower() for l in c_leaves) and
                                            not any(other_sp in l.name.lower() for l in c_leaves)
                                        )
                                        if only_target:
                                            return f"{c.dist:.6f}"
                                    # Descend dans les enfants qui ont encore target_sp
                                    for child in c.children:
                                        if any(target_sp in l.name.lower() for l in child.get_leaves()):
                                            queue_d.append(child)
                                return None

                            dist_sp1 = first_specific_dup(spec_node, sp1_name, sp2_name)
                            dist_sp2 = first_specific_dup(spec_node, sp2_name, sp1_name)

                            if dist_sp1 is not None and dist_sp2 is not None:
                                return dist_sp1, dist_sp2
                            return "N/A", "N/A"

###########BOUCLE

                        # On applique exactement les mêmes règles algorithmiques de 1 à 4 spéciations
                        for i in range(4):
                            if i < len(all_S_sorted):
                                current_spec = all_S_sorted[i]
                                
                                # 1) Duplications en AMONT (Pre-Spéciation i jusqu'à root_bloc)
                                nb_dup_pre = 0
                                curr_up = current_spec
                                while curr_up != root_bloc:
                                    curr_up = curr_up.up
                                    if curr_up and curr_up.name and curr_up.name.upper().startswith("D"):
                                        nb_dup_pre += 1
                                        
                                # 2) Duplications en AVAL (Post-Spéciation i - arrêt si autre S)
                                nb_dup_post = 0
                                queue = list(current_spec.children)
                                while queue:
                                    curr_down = queue.pop(0) 
                                    if not curr_down.is_leaf():
                                        if curr_down.name and curr_down.name.upper().startswith("S"):
                                            continue
                                        elif curr_down.name and curr_down.name.upper().startswith("D"):
                                            nb_dup_post += 1
                                            queue.extend(curr_down.children)
                                        else:
                                            queue.extend(curr_down.children)
                                            
                                pre_counts.append(nb_dup_pre)
                                post_counts.append(nb_dup_post)
                                #DETAILS
                                        # alors entre dans le bloc, regarde si c'est 1 S=> arret
                                        # si c'est un D => compte nb dup +1
                                        # si on continue de chercher on prend les enfants de ce neouds et on les ajoutent à la file d'attente "queue"
                                        # si curr_down.is_leaf est vrait le script lit "if not True" ce qui est faux=> il ignore le bloc et passe a la sequence suivante dans la file d'attente
                                        # pourquoi c'est OK :
                                        # car une feuille terminal ne sera jamais une duplication donc ne peut augmenter le compteur nb_dup_post
                                        # faire curr_down.children sur une feuille ferait planter le scrript ou list vide 
                                
                                # DISTANCE DES DUP 
                                # CORRECTION DU BUG down21: Appel explicite de la fonction et ajout des distances dans les listes
                                d1, d2 = get_robust_dup_distances(current_spec, "smaritimus", "salterniflorus")
                                SD_sp1_list.append(d1)
                                SD_sp2_list.append(d2)
                                    
                                # Grille de classification stricte
                                if nb_dup_pre >= 2:
                                    s_classifications.append("DFND")
                                elif nb_dup_pre < 2 and nb_dup_post >= 1:
                                    s_classifications.append("DFD")
                                elif nb_dup_pre ==0 and nb_dup_post == 0: #ajout 05/06 down 21.2
                                    s_classifications.append("absD")    
                                else:
                                    s_classifications.append("autres")
                            else:
                                # majdown21 Si la spéciation n'existe pas dans ce bloc, on remplit par N/A
                                pre_counts.append("N/A")
                                post_counts.append("N/A")
                                SD_sp1_list.append("N/A")
                                SD_sp2_list.append("N/A")
                                s_classifications.append("N/A")

################################################################################ #############
# CALCUL DE LA COLONNE SYNTHÉTIQUE GLOBALE "type_scenario"
################################################################################ #############
                        first2 = s_classifications[0:2]
                        has_DFND = "DFND" in first2
                        has_DFD = "DFD" in first2
                        has_absD = "absD" in first2 #ajout 11/06 down 21.2
                        
                        if has_DFND and has_DFD:
                            type_scenario = "mixte"
                        elif has_DFND:
                            type_scenario = "DFND"
                        elif has_DFD:
                            type_scenario = "DFD"
                        elif has_absD:
                            type_scenario = "absD" #ajout 05/06 down 21.2
                        else:
                            type_scenario = "autres"
                                
                        #############################################################################################
                        # SOUS TOTAL : SOMME DES ÉVÉNEMENTS PAR BLOC (LIGNE PAR LIGNE)
                        #############################################################################################
                        t_bloc_DFND = s_classifications.count("DFND")
                        t_bloc_DFD = s_classifications.count("DFD")
                        t_bloc_absD = s_classifications.count("absD") #ajout 05/06 down 21.2
                        t_bloc_autres = s_classifications.count("autres")        
                                
                        bloc_newick = root_bloc.write(format=1).replace('"', "'")
                       
##################################################################################################
# ÉTAPE 7 : ASSEMBLAGE CSV
##################################################################################################
                        row = [
                            f'"{og_id}"', f'"{bloc_id}"', f'"{smar_count}"', f'"{salter_count}"', f'"{seq_sum}"',
                            f'"{total_dup_bloc}"', f'"{total_spec_bloc}"', f'"{total_t_bloc}"', f'"{root_type}"', 
                            f'"{dist_1st_spec:.6f}"' if isinstance(dist_1st_spec, float) else f'"{dist_1st_spec}"', 
                            f'"{dist_2nd_spec:.6f}"' if isinstance(dist_2nd_spec, float) else f'"{dist_2nd_spec}"',
                            
                            f'"{pre_counts[0]}"', f'"{post_counts[0]}"', f'"{SD_sp1_list[0]}"', f'"{SD_sp2_list[0]}"',
                            f'"{pre_counts[1]}"', f'"{post_counts[1]}"', f'"{SD_sp1_list[1]}"', f'"{SD_sp2_list[1]}"',
                            f'"{pre_counts[2]}"', f'"{post_counts[2]}"', f'"{SD_sp1_list[2]}"', f'"{SD_sp2_list[2]}"',
                            f'"{pre_counts[3]}"', f'"{post_counts[3]}"', f'"{SD_sp1_list[3]}"', f'"{SD_sp2_list[3]}"',
                            
                            f'"{s_classifications[0]}"', f'"{s_classifications[1]}"', f'"{s_classifications[2]}"', f'"{s_classifications[3]}"', 
                            f'"{type_scenario}"',
                            f'"{t_bloc_DFND}"', f'"{t_bloc_DFD}"', f'"{t_bloc_absD}"', f'"{t_bloc_autres}"', 
                            f'"{bloc_newick}"'
                        ]
                        # Écriture physique de la ligne dans le fichier CSV, séparée par des virgules
                        f_out.write(",".join(row) + "\n")
                        # Incrémentation compteur de blocs écrits pour le rapport final
                        total_blocs_count += 1
                        # On mémorise qu'au moins un bloc a été extrait avec succès pour cet orthogroupe
                        og_written = True
                        # Passage à l'index de bloc suivant (ex: Bloc 2) si l'arbre cache d'autres bouquets distincts
                        bloc_id += 1
                # Si l'orthogroupe a généré au moins un bloc valide durant son parcours
                if og_written:
                    # On comptabilise cet OG dans notre indicateur de richesse des données
                    total_ogs_with_blocs += 1
            # Gestion des exceptions pour passer à l'arbre suivant si un fichier est totalement corrompu
            except Exception:
                # L'instruction continue permet d'ignorer le bug et de passer au fichier suivant de la boucle
                continue
##################################################################################################
# ÉTAPE 8 : STATISTIQUES FINALES GENERALE SLURM
##################################################################################################
        f_out.write("\n# === STATISTIQUES FINALES ===\n")
        f_out.write(f"# Nombre d'arbres de genes parcourus (fichiers lus) : {total_trees_scanned}\n")
        f_out.write(f"# Nombre d'OGs uniques presents dans le fichier (avec blocs) : {total_ogs_with_blocs}\n")
        f_out.write(f"# Nombre total de blocs (bouquets) ecrits : {total_blocs_count}\n")
        
# SECTION POUR LOGS 
        print("\n# === STATISTIQUES FINALES (Affichage Slurm) ===")
        print(f"Nombre d'arbres de gènes parcourus (fichiers lus) : {total_trees_scanned}")
        print(f"Nombre d'OGs uniques presents dans le fichier (avec blocs) : {total_ogs_with_blocs}")
        print(f"Nombre total de blocs (bouquets) ecrits : {total_blocs_count}\n")
        
##################################################################################################
# ÉTAPE 9 : SÉPARATION EN 4 SOUS-CSV BASÉE SUR type_scenario & CALCUL DU TABLEAU RÉSUMÉ
##################################################################################################    
       
    # Déterminer => noms des 4 fichiers de sortie basés sur le nom du CSV global
    base_name, ext = os.path.splitext(output_csv)
    csv_DFND = f"{base_name}_DFND{ext}"
    csv_DFD = f"{base_name}_DFD{ext}"
    csv_mixte = f"{base_name}_mixte{ext}"
    csv_absD = f"{base_name}_absD{ext}"
    csv_autres = f"{base_name}_autres{ext}"
    
    # Initialisation du compteur @ résumé
    summary_counts = {"DFND": 0, "DFD": 0, "mixte": 0, "absD": 0, "autres": 0}

    # Ouverture des fichiers en écriture
    with open(output_csv, "r") as f_in, \
         open(csv_DFND, "w") as f_DFND, \
         open(csv_DFD, "w") as f_DFD, \
         open(csv_mixte, "w") as f_mix, \
         open(csv_absD, "w") as f_absD, \
         open(csv_autres, "w") as f_aut:
         
        # Récupération de la ligne d'en-tête et écriture dans les 4 fichiers
        header_line = f_in.readline()
        f_DFND.write(header_line)
        f_DFD.write(header_line)
        f_mix.write(header_line)
        f_absD.write(header_line)
        f_aut.write(header_line)
        
       # Trouver l'index de 'type_scenario' dynamiquement SPETdown21
        header_cols = header_line.strip().split(",")
        try:
            idx_scenario = header_cols.index("type_scenario")
        except ValueError:
            print("ERREUR : Colonne 'type_scenario' introuvable dans le fichier CSV.")
            return
        
        # Parcours de chaque ligne du CSV général
        for line in f_in:
            # Si on atteint le bloc des statistiques finales, on arrête la lecture
            if line.startswith(("#", "\n")):
                continue
                
            # Séparation des colonnes par la virgule
            parts = line.strip().split(",")
            # CORRECTION INDEX DYANMIQUE SPETdown21
            if len(parts) > idx_scenario:
                classif1 = parts[idx_scenario].replace('"', '').strip()
                
                # Distributions
                if classif1 == "DFND":
                    f_DFND.write(line)
                    summary_counts["DFND"] += 1
                elif classif1 == "DFD":
                    f_DFD.write(line)
                    summary_counts["DFD"] += 1
                elif classif1 == "mixte":
                    f_mix.write(line)
                    summary_counts["mixte"] += 1
                elif classif1 == "absD":
                    f_absD.write(line)
                    summary_counts["absD"] += 1
                else:
                    f_aut.write(line)
                    summary_counts["autres"] += 1
################################################################################################
# ÉTAPE 10 : AFFICHAGE RÉSUMÉ 
###########################################################################################
    with open(output_csv, "a") as f_out_append:
        f_out_append.write("\n# === TABLEAU RÉSUMÉ DES SCÉNARIOS PAR BLOC ===\n")
        f_out_append.write(f"# Catégorie DFND  : {summary_counts['DFND']} bloc(s)\n")
        f_out_append.write(f"# Catégorie DFD    : {summary_counts['DFD']} bloc(s)\n")
        f_out_append.write(f"# Catégorie mixte : {summary_counts['mixte']} bloc(s)\n")
        f_out_append.write(f"# Catégorie absD : {summary_counts['absD']} bloc(s)\n")
        f_out_append.write(f"# Catégorie autres: {summary_counts['autres']} bloc(s)\n")

    print("# === TABLEAU RÉSUMÉ DES SCÉNARIOS PAR BLOC (Affichage Slurm) ===")
    print(f"Nombre de blocs DFND   : {summary_counts['DFND']}")
    print(f"Nombre de blocs DFD     : {summary_counts['DFD']}")
    print(f"Nombre de blocs mixtes : {summary_counts['mixte']}")
    print(f"Nombre de blocs absD : {summary_counts['absD']}")
    print(f"Nombre de blocs autres : {summary_counts['autres']}\n")

    print("# === SÉPARATION DES COMPATIBILITÉS TERMINÉE ===")
    print(f"Fichier DFND créé : {csv_DFND}")
    print(f"Fichier DFD créé  : {csv_DFD}")
    print(f"Fichier mixte créé        : {csv_mixte}")
    print(f"Fichier absD créé        : {csv_absD}")
    print(f"Fichier autres créé       : {csv_autres}\n")

# Condition standard en Python pour lancer la fonction main() uniquement si le script est exécuté directement
if __name__ == "__main__":
    main()
