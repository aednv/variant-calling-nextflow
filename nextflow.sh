#!/bin/bash

#BSUB -n 1
#BSUB -J nextflow_job
#BSUB -o nextflow_job.out
#BSUB -e nextflow_job.err
#BSUB -q long
#BSUB -W 30:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"

module load nextflow/20.10.0.5430

nextflow run call-variants.nf -with-dag results/flowchart.html -resume
 
