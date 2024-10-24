#!/bin/bash
#SBATCH --partition=short
#SBATCH --ntasks=10
#SBATCH --mem=60G
#SBATCH --time=0-06:00:00

module load spaceranger/2.1.0


cd /project/shared/spatial_data_camp/spaceranger/outputs

fastq=/project/shared/spatial_data_camp/spaceranger/fastq
reference=/project/shared/spatial_data_camp/spaceranger/reference/refdata-gex-GRCh38-2020-A
image=/project/shared/spatial_data_camp/spaceranger/images/A1_colon_day0.tif


spaceranger count --id=SRR14083626 \
                   --transcriptome=$reference \
                   --fastqs=$fastq \
                   --sample=SRR14083626 \
                   --image=$image \
                   --slide V19S23-097 \
                   --area=A1 \
                   --localcores=10 \
                   --localmem=60 \
                   --reorient-images true

