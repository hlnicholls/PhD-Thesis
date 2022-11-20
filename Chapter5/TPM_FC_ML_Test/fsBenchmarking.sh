#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 4      # Request 1 core
#$ -l h_rt=72:0:0 # Request 24 hour runtime
#$ -l h_vmem=20G   # Request 1GB RAM

module load python

source /data/WHRI-Bioinformatics/Nicholls/python/BPGWASPredict/BPGWASPredict_env/bin/activate

python FSML_pipeline.py