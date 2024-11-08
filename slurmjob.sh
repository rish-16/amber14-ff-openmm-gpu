#!/bin/bash

#SBATCH --job-name=grassy_12k_processed
#SBATCH --time 10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=gpu
#SBATCH --gpus=2
#SBATCH --mem=40G
#SBATCH --output=./logs/slurm/%x_%j.out
#SBATCH --error=./logs/slurm/%x_%j.err
cd /gpfs/gibbs/pi/krishnaswamy_smita/ra852/GRASSY
module load miniconda
conda activate grassy

python src/train.py trainer=ddp
# salloc --gpus=2 --time=2:00:00 --partition gpu_devel