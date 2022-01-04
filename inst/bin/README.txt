Script to make test VCFs using ANGIOPREDICT data, contact bruce01campus@gmail.com

Command:

nextflow run main.nf -c ~/.nextflow/assets/brucemoran/batch_somatic/nextflow.config -profile genome,singularity --refDir /data/genome/bmoran/reference/somatic_n-of-1 --bedInp "/store2/bmoran/ANGIOPREDICT/exome_APD_rerun/BRAF_KRAS/BRAF_KRAS.bed" --runID "BRAF_KRAS_test"
