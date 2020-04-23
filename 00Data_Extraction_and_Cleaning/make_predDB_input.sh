#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb

zcat /home/ryan/topmed/multiomic_modeling/predictor_matrices/dosages/topmed_dosages/chr${chr}.maf0.01.R20.8.dosage.txt.gz | cut -f 2,7- -d " " >> /home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr${chr}_geno.txt
zcat /home/ryan/topmed/multiomic_modeling/predictor_matrices/dosages/topmed_dosages/chr${chr}.maf0.01.R20.8.dosage.txt.gz | awk '{varid=$1 "_" $3 "_" $4 "_" $5 "_b37" ;print $1 " " $3 " " varid " " $3 " " $4 " " $5 " " $2}' >> /home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr${chr}_anno.txt
sed -i '1 s/[0-9]\+_//g' /home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr${chr}_geno.txt
sed -i '1 s/chr_pos_ref_alt_b37/varid/g' /home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr${chr}_anno.txt
gzip -f /home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr${chr}_geno.txt
gzip -f /home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr${chr}_anno.txt
