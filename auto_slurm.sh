#!/bin/bash
#SBATCH --job-name=hyper-vlasov
#SBATCH --output=log-%j.txt
#SBATCH --error=errors-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=16G
#SBATCH --time=2:00:00

srun bin/hyper-vlasov output/out.vtk