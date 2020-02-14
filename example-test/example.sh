#!/bin/bash

# LOHHLA_dir="~/libs/LOHHLA"
LOHHLA_dir='/home/m.slagter/libs/LOHHLA'
novo_loc='/home/m.slagter/conda/bin'
gatk_loc='/DATA/users/l.fanchi/libs/picard-tools-1.120'
jellyfish_loc='/home/m.slagter/conda/bin/jellyfish'

Rscript "$LOHHLA_dir/LOHHLAscript.R" \
  --patientId example \
  --outputDir $LOHHLA_dir/example-test/out/example-out/ \
  --normalBAMfile $LOHHLA_dir/example-test/bam/example_BS_GL_sorted.bam \
  --tumorBAMfile $LOHHLA_dir/example-test/bam/example_tumor_sorted.bam \
  --BAMDir $LOHHLA_dir/example-test/bam/  \
  --hlaPath $LOHHLA_dir/example-test/hlas \
  --HLAfastaLoc $LOHHLA_dir/data/example.patient.hlaFasta.fa \
  --HLAexonLoc $LOHHLA_dir/data/hla.dat \
  --CopyNumLoc $LOHHLA_dir/example-test/solutions.txt \
  --genomeAssembly hg19 \
  --mappingStep TRUE \
  --kmerSize 50 \
  --plottingStep FALSE \
  --minCoverageFilter 10 \
  --forceRedo TRUE \
  --fishingStep TRUE \
  --cleanUp FALSE \
  --gatkDir $gatk_loc \
  --novoDir $novo_loc
