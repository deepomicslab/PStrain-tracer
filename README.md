# PStrain-tracer
## What is PStrain-tracer?
PStrain-tracer is a downstream analysis and visulization package of [**PStrain**](https://github.com/wshuai294/PStrain), a tool to infer the genome sequence and propotion of strains from whole genome shotgun metagenomics data.
## What could PStrain-tracer do?
![flowchart](PStrain-tracer.png)
## Installation
PStrain-tracer was written in python3 and perl5, and using R code to do the visualiztion.
```
  git clone https://github.com/deepomicslab/PStrain-tracer.git
```
- Python Dependencies
  - numpy
- Configuration
  - modify the ```config.list``` in the ```src/``` folder
```
Rscript=/where/R/installed/bin/Rscript
R_lib=/where/R/library

MUMmer=/where/MUMmer/install/MUMmer3.23/
fasttree=/where/excutable/FastTree
ass_sum=/where/PStrain-tracer/src/assembly_summary_refseq_20200423.txt
species_marker=/where/PStrain/db/species_markers.txt.gz
marker_gene=/where/PStrain/db/marker_gene.fna.gz
```
If you don't use the module of report the most similar strains from Genbank, you only need modify Rscript and R_lib.
- The R_lib should contains:
  - ggplot2
  - RColorBrewer
  - pheatmap
- For the report the similar strains modules, you need:
  - [**MUMmer**](https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download)
  - [**FastTree2**](http://www.microbesonline.org/fasttree/)
  - [**Assembly Summary**](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt) from NCBI Assemebly FTP.
  - or you can use the ```assembly_summary_refseq_20200423.txt``` we provided in ```src``` folder, but remember gunzip before use.
  - ```species marker``` and ```marker gene``` are provided in [**PStrain**](https://github.com/wshuai294/PStrain).

## Example
We provide an example to run the package in ```example``` folder. ```example.sh``` showed the commands used in a compeleted work.
After your installation and configuration. You can change your direction to the ```example``` folder and directly run ```example.sh``` to figure out if it worked well.
```
# decompressed the demo PStrain inputs
cat demo_PStrain_result.tar.gz.* | tar -zxv 

# Cluster strains based on genotypes
python ../src/PT-01-cluster.py -c config -o demo_PStrain_result --similarity 0.9

# Create the output folder
mkdir output

perl ../src/PT-02-FMT_visual.pl  -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9

perl ../src/PT-03-engraft.pl -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9
perl ../src/PT-03-engraft.pl -WDR output -LS experimets_list.tsv -PS demo_PStrain_result/ -S 0.9 -m abd

perl ../src/PT-04-heterzygousSNV_strainN.pl -WDR output -PS demo_PStrain_result/ -S 0.9
perl ../src/PT-05-wilcox_testing.pl -WDR output -PS demo_PStrain_result/ -S 0.9 -LS group_info.tsv
perl ../src/PT-06-detect.pl -WDR output -PS demo_PStrain_result/ -S 0.9 -LS strain.txt
```

## Output

