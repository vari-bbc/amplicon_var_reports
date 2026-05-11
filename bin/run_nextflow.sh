#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J nextflow_workflow
#SBATCH -o nextflow_workflow.o
#SBATCH -e nextflow_workflow.e
#SBATCH --ntasks 1
#SBATCH --time 120:00:00
#SBATCH --mem=8G
#SBATCH --partition=long


nextflow_module="bbc2/nextflow/nextflow-25.04.07" 

module load $nextflow_module

nextflow run main.nf -profile vai
