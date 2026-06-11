#!/bin/bash
#SBATCH --job-name=nb_bouquet
#SBATCH -o /home/genouest/cnrs_umr6553/aubonnier/EXP/output/%x.%j.out
#SBATCH -e /home/genouest/cnrs_umr6553/aubonnier/EXP/log/%x.%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 

echo -e "\t start script (step1) \n"

# Sécurité Bash
#set -euo pipefail

echo -e "\t load env python (step2) \n"

# Chargement de l'environnement Python Genouest
. /local/env/envpython-3.9.5.sh
echo -e "\t test version python (step3) \n"
python --version
echo -e "\t preparatifs termines (step4) \n"

# Chargement chemins
if [ -f "/home/genouest/cnrs_umr6553/aubonnier/scriptab/pathAB.sh" ]; then
    source /home/genouest/cnrs_umr6553/aubonnier/scriptab/pathAB.sh
else
    echo "WARNING: pathAB.sh non trouvé"
    exit 1
fi

echo -e "\t Variables chargées (step5) \n"

# Vérif variable critique
[ -z "${scriptab:-}" ] && { echo "ERREUR: scriptab non défini"; exit 1; }

echo -e "\t Test variables critique (step6) \n"

outfile="/home/genouest/cnrs_umr6553/aubonnier/EXP/output/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.csv"

echo -e "\t Outfile par défaut définir (step7) \n"

# Récupération des arguments
input=""

echo -e "\t Chargement de l'argument input (step8) \n"

while getopts "i:o:" opt; do
    case $opt in
        i) input="$OPTARG" ;;
        o) outfile="$OPTARG" ;;  # écrase le défaut si fourni
        *) exit 1 ;;
    esac
done

echo -e "\t Lecture des arguments par getopts (step9) \n"

# Vérification des variables obligatoires
[ -z "${input}" ] && { echo "ERREUR: Dossier input manquant (-i)"; exit 1; }
[ -z "${outfile}" ] && { echo "ERREUR: Fichier output manquant (-o)"; exit 1; }
echo -e "\t Vérification des variables obligatoire validée (step10) \n"

echo "Start: $(date)"
echo -e "\t Date démarrage script annexe (step11) \n"

mkdir -p "$(dirname "$outfile")"
echo -e "\t Création dossier outfile (step12) \n"

# Retrait inutile Montrer fichier de résultats retiré
# echo "OG,Smaritimus_count,Salterniflorus_count" > "$outfile"

#Activation conda
. /local/env/envconda.sh
echo -e "\t Chargement outils conda genouest(step13) \n"
eval "$(conda shell.bash hook)"
conda activate ete3_env
echo -e "\t Chargement envab ete3 validé(step14) \n"

python3 -c "import ete3; print('ETE3 OK')"

# appel au script python pour compter via ETE3
python3 "$scriptab/BM2_tempo2.py" "$input" "$outfile"

conda deactivate
echo -e "\t Conda désactivé (step15) \n"

echo -e "\t Fin script annexe (step16) \n"

echo "Finish: $(date)"
echo "Résultats enregistrés dans : $outfile"
