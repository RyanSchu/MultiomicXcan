#!/bin/bash
#PBS -S /bin/bash
#PBS -o /home/ryan/software/TOPMED/Proteome/impute_expression_analysis/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/software/TOPMED/Proteome/impute_expression_analysis/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb
p=$1
/usr/local/bin/MetaXcan_software/PrediXcan.py \
--throw \
--text_genotypes /home/ryan/topmed/proteome/${p}/genotypes/unfiltered_dosages/topmed_dosages/predixcan_input/chr*.maf0.0.R20.8.dosage.txt.gz \
--text_sample_ids /home/ryan/topmed/proteome/${p}/genotypes/unfiltered_dosages/topmed_dosages/prediXcan_samples.txt \
--model_db_path /home/ashley/TopMed/EN_Base_Models/dbs/CAU_WG_unfiltered.db \
--prediction_output /home/ryan/topmed/multiomic_modeling/imputed_results/CAU_baseline_model_${p}_geno_predicted_expression.txt \
--prediction_summary_output /home/ryan/topmed/multiomic_modeling/imputed_results/CAU_baseline_model_${p}_geno_prediction_summary.txt

echo /home/ryan/topmed/multiomic_modeling/imputed_results/CAU_multiomic_model_${p}_geno_prediction_summary.txt
