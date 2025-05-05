#!/bin/bash
#SBATCH --mem=3000
#SBATCH --time=24:00:00 --qos=1day #adjust when increasing sample size
#SBATCH --job-name=gBLUISA
#SBATCH --cpus-per-task=1   #make sure to modify $cpus too!!!
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --array=1-2000

#8245 genes

set -e
date

#where are the files
gemmafiles=/Genomics/ayroleslab2/lamaya/bigProject/Genotypes_feb2020/headbody_notsameindv/ForGemma

#make fam into a file with expressio info - but do this before running the array, otherwise 
#may jobs trying to do the same at the same ime generates errors

grm=${gemmafiles}/output/ctrl_headbody.final.bodyonlynorelatedness1.LUISAind.Sep221.miss50.sXX.txt

#file with the order of genes in the fam file
geneorder=${gemmafiles}/geneorderFAMfile_body_LUISAind.Sep221.txt
	#geneorder=${gemmafiles}/geneorderFAMfile_199genes.txt
#gene id forthe output files
#the max array number is 2000, so I have to set up an offset to be able to read all genes

#covariate file with first column being 1 representing the intersect
cov=${gemmafiles}/ctrl_headbody.final.body.covs.forGemma.LUISAind.Sep221

#plink files
bfile=ctrl_headbody.final.bodyonlynorelatedness1.LUISAind.Sep221


###Round 1

offset=0    #2000 #4000 #6000 #8000 (when doing 8000 set array from 1-277)
gene=`awk -v file=$SLURM_ARRAY_TASK_ID -v off=$offset '{if (NR==file+off) print $0 }' ${geneorder}`

echo "${gene} being mapped" 
echo "offset ${offset}"
geneline=$(($offset+$SLURM_ARRAY_TASK_ID))

~/bin/gemma-0.98.1-linux-static -bfile ${gemmafiles}/${bfile} -k ${grm} -maf 0.05 -miss 0.5 -lmm 1 -n $geneline -o indLUISA/${gene}.LUISAind -c ${cov}

#make table smaller to load into R 

awk '{ if ( $12 < 0.01 ) print $0 }' output/indLUISA/${gene}.LUISAind.assoc.txt > output/indLUISA/${gene}.LUISAind.assoc.small.txt
sed -i '1ichr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald'  output/indLUISA/${gene}.LUISAind.assoc.small.txt

echo 'Done!'

### Round 2

offset=2000    #2000 #4000 #6000 #8000 (when doing 8000 set array from 1-277)
gene=`awk -v file=$SLURM_ARRAY_TASK_ID -v off=$offset '{if (NR==file+off) print $0 }' ${geneorder}`

echo "${gene} being mapped"
echo "offset ${offset}"
geneline=$(($offset+$SLURM_ARRAY_TASK_ID))

~/bin/gemma-0.98.1-linux-static -bfile ${gemmafiles}/${bfile} -k ${grm} -maf 0.05 -miss 0.5 -lmm 1 -n $geneline -o indLUISA/${gene}.LUISAind -c ${cov}

#make table smaller to load into R

awk '{ if ( $12 < 0.01 ) print $0 }' output/indLUISA/${gene}.LUISAind.assoc.txt > output/indLUISA/${gene}.LUISAind.assoc.small.txt
sed -i '1ichr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald'  output/indLUISA/${gene}.LUISAind.assoc.small.txt

echo 'Done!'

### Round 3

offset=4000    #2000 #4000 #6000 #8000 (when doing 8000 set array from 1-277)
gene=`awk -v file=$SLURM_ARRAY_TASK_ID -v off=$offset '{if (NR==file+off) print $0 }' ${geneorder}`

echo "${gene} being mapped"
echo "offset ${offset}"
geneline=$(($offset+$SLURM_ARRAY_TASK_ID))

~/bin/gemma-0.98.1-linux-static -bfile ${gemmafiles}/${bfile} -k ${grm} -maf 0.05 -miss 0.5 -lmm 1 -n $geneline -o indLUISA/${gene}.LUISAind -c ${cov}

#make table smaller to load into R

awk '{ if ( $12 < 0.01 ) print $0 }' output/indLUISA/${gene}.LUISAind.assoc.txt > output/indLUISA/${gene}.LUISAind.assoc.small.txt
sed -i '1ichr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald'  output/indLUISA/${gene}.LUISAind.assoc.small.txt

echo 'Done!'

### Round 4

offset=6000    #2000 #4000 #6000 #8000 (when doing 8000 set array from 1-277)
gene=`awk -v file=$SLURM_ARRAY_TASK_ID -v off=$offset '{if (NR==file+off) print $0 }' ${geneorder}`

echo "${gene} being mapped"
echo "offset ${offset}"
geneline=$(($offset+$SLURM_ARRAY_TASK_ID))

~/bin/gemma-0.98.1-linux-static -bfile ${gemmafiles}/${bfile} -k ${grm} -maf 0.05 -miss 0.5 -lmm 1 -n $geneline -o indLUISA/${gene}.LUISAind -c ${cov}

#make table smaller to load into R

awk '{ if ( $12 < 0.01 ) print $0 }' output/indLUISA/${gene}.LUISAind.assoc.txt > output/indLUISA/${gene}.LUISAind.assoc.small.txt
sed -i '1ichr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald'  output/indLUISA/${gene}.LUISAind.assoc.small.txt

echo 'Done!'

### Round 5

offset=8000    #2000 #4000 #6000 #8000 (when doing 8000 set array from 1-277)
gene=`awk -v file=$SLURM_ARRAY_TASK_ID -v off=$offset '{if (NR==file+off) print $0 }' ${geneorder}`

echo "${gene} being mapped"
echo "offset ${offset}"
geneline=$(($offset+$SLURM_ARRAY_TASK_ID))

~/bin/gemma-0.98.1-linux-static -bfile ${gemmafiles}/${bfile} -k ${grm} -maf 0.05 -miss 0.5 -lmm 1 -n $geneline -o indLUISA/${gene}.LUISAind -c ${cov}

#make table smaller to load into R

awk '{ if ( $12 < 0.01 ) print $0 }' output/indLUISA/${gene}.LUISAind.assoc.txt > output/indLUISA/${gene}.LUISAind.assoc.small.txt
sed -i '1ichr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald'  output/indLUISA/${gene}.LUISAind.assoc.small.txt


date
echo 'Done!'
