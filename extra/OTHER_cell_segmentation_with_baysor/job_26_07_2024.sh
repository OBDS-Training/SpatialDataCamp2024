#!/bin/bash
#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=20
#SBATCH --mem=900G
#SBATCH --partition=long


export JULIA_NUM_THREADS=20 


out=/project/shared/spatial_data_camp/datasets/PRECOMPUTED/baysor
transcripts=/project/shared/spatial_data_camp/notebooks_and_code/EXTRA/OTHER_cell_segmentation_with_baysor/detected_transcripts_for_baysor.csv

baysor run $transcripts :cell_id \
                        -m 10 \
                        -s 6.5 \
                        -x x_location \
                        -y y_location \
                        -z z_location \
                        -g feature_name \
                        --save-polygons=geojson \
                        -p \
                        -o $out \
                        --n-clusters 12 \
                        --prior-segmentation-confidence 0.5 \
                        --count-matrix-format tsv

