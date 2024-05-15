#!/bin/bash
#SBATCH --job-name=2D_intm1        # nom du job
#SBATCH --time=2-00:00:00          # temps maximal d'exécution du job
#SBATCH --mail-user=victor.ballester-ribo@dauphine.eu --mail-type=ALL
#SBATCH --clusters=mesopsl1
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --account=vballester
#SBATCH --nodes=3                # Nombre total de processus MPI
#SBATCH --ntasks-per-node=16       # Nombre de processus MPI par noeud
#SBATCH --hint=nomultithread       # 1 processus MPI par coeur physique (pas d'hyperthreading)
#SBATCH --output=TravailMPI%j.out  # Nom du fichier de sortie
#SBATCH --error=TravailMPI%j.out   # Nom du fichier d'erreur (ici commun avec la sortie)

# on se place dans le répertoire de soumission
cd ${SLURM_SUBMIT_DIR}

# nettoyage des modules charges en interactif et herites par defaut
module purge

# chargement des modules
module load intel-all/19.0.4

# echo des commandes lancées
set -x

# exécution du code
srun ./hd2D
