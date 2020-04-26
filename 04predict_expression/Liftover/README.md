For this example we will use the MESA dbs from predictdb.org 

Scripts and tools for lifting over db file/prediXcan weights file. Bear in mind that doing lifting over models between builds may not be equivalent to building models in each build separately. You should get most of the same snps, but between builds snps change position and asome get deleted, split, or categorized differently. So it could be that in one build a snp within the cis window for a model, but in another build it falls out of the cis window. It could be that one build has removed a snp contained by another, or what was previously considered an indel has been recategorized as several single point mutations. For best results you should really build models in each build that you're considering, but for the most part it is adequate to liftover the snps in one build or another.

## SNP ids
Performing liftover on models is only truly necessary when the models have ids in cpos format. If all files are in rsid format, this process is redundant as the models can map to the correct snp based on the rsid. In cpos however, these snps will have different ids depending on the build and cannot map to each other between builds.

first extract the weights using the `sql.py` script. 

```
python sql.py --db ALL_imputed_10_peer_3_pcs_v2.db --out ALL_imputed_10_peer_3_pcs_v2_weights.txt
```

Then put the weights into liftover format

```
awk '{split($3,a,"_"); print "chr" a[1] FS a[2] FS a[2] + 1 FS $3}' ALL_imputed_10_peer_3_pcs_v2_weights.txt > prelift_ALL.txt
```

Then liftover

```
 liftOver prelift_ALL.txt /home/rschubert1/data/liftover/hg19ToHg38.over.chain.gz lifted_hg38_ALL.txt unmapped_ALL.txt
```

Next generate the cpos ids and reformat the varids, also dropping alt and random chromosomes

```
awk '{split($4,a,"_");gsub("chr","",$1); print $0 FS "hg38:" $1 ":" $2 FS $1 "_" $2 "_" a[3] "_" a[4] "_b38"}' lifted_hg38_ALL.txt |grep -v -F alt | grep -v -F rand > lifted_hg38_ALL_with_Varid.txt
```

Then use the `liftover_weights.R` script to make a new weights table with updated ids

```
Rscript liftover_weights.R --lifted lifted_hg38_ALL_with_Varid.txt --weights ALL_imputed_10_peer_3_pcs_v2_weights.txt --out lifted_ALL_imputed_10_peer_3_pcs_v2_weights.txt
```

now finally you can run `ALL_sql_eof` to update you table

```
cp ALL_imputed_10_peer_3_pcs_v2.db lifted_ALL_imputed_10_peer_3_pcs_v2.db
bash ALL_sql_eof
```
