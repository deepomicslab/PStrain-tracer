# decompressed the demo PStrain inputs
cat demo_PStrain_result.tar.gz.* | tar -zxv 

python ../src/PT-01-cluster.py -c config -o demo_PStrain_result --similarity 0.9

mkdir output

perl ../src/PT-02-FMT_visual.pl  -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9

perl ../src/PT-03-engraft.pl -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9 -m relative
perl ../src/PT-03-engraft.pl -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9

perl ../src/PT-04-heterzygousSNV_strainN.pl -WDR output  -PS demo_PStrain_result/ -S 0.9
perl ../src/PT-05-wilcox_testing.pl -WDR output -PS demo_PStrain_result/ -S 0.9 -LS group_info.tsv
perl ../src/PT-06-detect.pl -WDR output -PS demo_PStrain_result/ -S 0.9 -LS strain.txt
