#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb

bcftools view \
/home/ryan/topmed/cau_super_pop/laurens_imputation/chr${chr}.dose.vcf.gz \
--samples-file /home/ryan/topmed/multiomic_modeling/TRUE_intersect_CAU_PBMC_FID_IID_list.txt \
-o /home/ryan/topmed/multiomic_modeling/predictor_matrices/vcfs/CAU_PBMC_chr${chr}.dose.vcf.gz \
-O z
