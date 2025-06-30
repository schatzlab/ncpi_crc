#!/bin/bash

### USAGE: ./run_ohana.sh <VCF> <IDENT> <DOWNSAMP_PATH>
### <VCF>: path to VCF of variants
### <IDENT>: VCF ID (ex: chromosome number)
### <DOWNSAMP_PATH>: path to dir with downsampled matrices 
###
### Runs Ohana on data from 1000 Genomes Phase 3
### using Q/C matrices inferred from k=8, mi=50

OHANAPATH=~/progs/ohana/bin # path to Ohana bin folder

VCF=$1
IDENT=$2
DOWNSAMP_PATH=$3

### Scan for frquency differentiation outliers:

# convert VCF into .ped file
plink --vcf ${VCF} \
      --out vcfs/${IDENT} \
      --recode 12 \
      --geno 0.0 \
      --tab

# convert .ped file into .dgm file (a.k.a. G matrix) for qpas
${OHANAPATH}/convert ped2dgm vcfs/${IDENT}.ped vcfs/${IDENT}.dgm
echo "vcf converted into dgm for "${IDENT}

# produce admixture-corrected AFs for full dataset by reading in pre-calculated Q matrix
# -fq: optimizer should not optimize the specific Q matrix (because it was done already, I assume)
${OHANAPATH}/qpas vcfs/${IDENT}.dgm \
                                            -k 8 \
                                            -qi ${DOWNSAMP_PATH}/chr21_pruned_50_Q.matrix \
                                            -fo f_matrices/${IDENT}_F.matrix \
                                            -e 0.0001 \
                                            -fq \
                                            -mi 50
echo "Finished inferring ancestral AFs for "${IDENT}

# apply selscan
for i in {1..8}
do
	matrix_num=$i
	${OHANAPATH}/selscan vcfs/${IDENT}.dgm \
		f_matrices/${IDENT}_F.matrix \
		${DOWNSAMP_PATH}/c_matrices/chr21_pruned_50_C.matrix \
		-cs ${DOWNSAMP_PATH}/c_matrices/chr21_pruned_50_C_p${matrix_num}.matrix \
		> selscan/${IDENT}_ohanascan_k8_p${i}.out
	echo "Finished running selscan on ancestry component "$i
done
