# decompressed the demo PStrain inputs
cat demo_PStrain_result.tar.gz.* | tar -zxv 

# Cluster strains based on genotypes
python ../src/PT-01-cluster.py -c config -o demo_PStrain_result --similarity 0.9

# Create the output folder
mkdir output
# Strain engraftment pattern and visualization
perl ../src/PT-02-FMT_visual.pl  -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9

# Engrafment summary and visualization, default is run with proportion, with '-m abd' will run with original input abundance
perl ../src/PT-03-engraft.pl -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9
perl ../src/PT-03-engraft.pl -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9 -m abd

# Number of strains and heterozygous SNVs statistic
perl ../src/PT-04-heterzygousSNV_strainN.pl -WDR output -PS demo_PStrain_result/ -S 0.9

# Detect potential determinant strains with wilcox testing
perl ../src/PT-05-wilcox_testing.pl -WDR output -PS demo_PStrain_result/ -S 0.9 -LS group_info.tsv

# Reporting the most similar genome of the provided strains
perl ../src/PT-06-detect.pl -WDR output -PS demo_PStrain_result/ -S 0.9 -LS strain.txt
