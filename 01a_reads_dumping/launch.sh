#!/bin/bash -l
 
################### Gaussian Job Batch Script Example ###################
# Section for defining queue-system variables:
#-------------------------------------
# SLURM-section
#SBATCH --job-name=Genomics
#SBATCH --nodes=1
#SBATCH --partition=all
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=sal_%j.out
#SBATCH --error=sal_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aaron.garcia@uclm.es
######################################
# Section for defining job variables and settings:


module load sratoolkit conda
 
export SCRDIR=/scratch/$USER/$SLURM_JOB_ID
 
mkdir -p "$SCRDIR"
 
submitdir=$SLURM_SUBMIT_DIR
tempdir=$SCRDIR
input="rnaseq"

# copia los ficheros de datos en $tempdir 
cp -r "$submitdir"/$input "$tempdir"
cd "$tempdir"/"$input"
# ejecuta tus comandos

conda activate rnaseq

while IFS= read -r sample_id; do
	bash pipeline.sh "$sample_id"
done < SRR_Acc_List.txt

conda deactivate

# recupera ficheros
cp -r "$SCRDIR"/"$input" "$submitdir"
 
# limpieza 
cd $submitdir
rm -r "$tempdir"/*
rmdir  "$tempdir"
 
echo "Job finished at"
date
################### Job Ended ###################
exit 0

