#!/usr/bin/perl -w
use strict;

use IO::Uncompress::Bunzip2 '$Bunzip2Error';

use Cwd qw/abs_path/;
use Getopt::Long;
use FindBin qw/$RealBin $RealScript/;

my ($whole_dir,$spe,$metaphlan_v,$threshold);
my ($Help);
my %DB;
my @software = (
    "fasttree", "nucmer","delta-filter","show-snps"
);
$threshold = 50;

print "Script directory: $RealBin\n";
my $database_list = "$RealBin/config.list";
$DB{'SGB2GTDB'} = "$RealBin/mpa_vOct22_CHOCOPhlAnSGB_202212_SGB2GTDB.tsv";
$DB{'gtdb'} = "$RealBin/bac120_metadata_r207.tsv";


GetOptions(
	# required para
	"WDR=s"	=>	\$whole_dir,
	"S=s"	=>	\$spe,
	"V=s"  =>	\$metaphlan_v,
	"I=s"	=>	\$DB{'in'},
    "N=s"   =>  \$threshold,
    "DBS=s"   =>  \$DB{'species_marker'},
    "DBM=s"   =>  \$DB{'marker_gene'},
    #"DBG=s"   =>  \$DB{'SGB2GTDB'},
	# help
	"h"		=>	\$Help
);


my $tmp_dir = "$whole_dir/find_strain/$spe";
print "Working directory: $tmp_dir\n";
# check para
&para_alert;


`zcat $DB{'species_marker'}|grep -w $spe > $tmp_dir/$spe.list`;
open IN,"$tmp_dir/$spe.list";

my %in = ();
while(<IN>){
    chomp;
    my @l = split /\s/;
    $in{$l[0]} = 1;   
}
open OU,">$tmp_dir/$spe.marker.fa";

my $in;
if ($DB{'marker_gene'} =~ /gz$/) {
  open my $fh, "-|", "gunzip -dc $DB{'marker_gene'}";
  $in = $fh;
} elsif ($DB{'marker_gene'} =~ /bz2$/) {
  $in = IO::Uncompress::Bunzip2->new($DB{'marker_gene'});  
} elsif ($DB{'marker_gene'} =~ /fna$/) {
  open my $fh, "<", $DB{'marker_gene'};
  $in = $fh;
}

$/ = '>';

<$in>;
while(my $seq=<$in>){

    $seq =~ s/>$//;
    my $id = (split /\s/,$seq)[0];
    if(exists $in{$id}){
        print OU ">$seq";
    }
}
$/ = "\n";

open SH1,">$tmp_dir/1.dl.sh";
`mkdir -p $tmp_dir/dl/`;

open SH2,">$tmp_dir/2.snp.sh";
`mkdir -p $tmp_dir/snp/`;

my %search_id =();
$search_id{'M3'} = 'ncbi_taxid';
$search_id{'M4'} = 'gtdb_taxonomy';


my $taxid = '';
if($metaphlan_v eq 'M3'){
    $taxid = `cat $tmp_dir/$spe.list|awk -F '__' '{print \$1}'|head -1`;
    chomp $taxid;
}elsif($metaphlan_v eq 'M4'){
    my $SGBid = `cat  $tmp_dir/$spe.list|awk '{print \$1}'|awk -F '|' '{print \$NF}'|head -1 `;
    chomp $SGBid;
    $taxid = `cat $DB{'SGB2GTDB'}|grep -w $SGBid|awk -F '\t' '{print \$2}'`;
    chomp $taxid;
}
if($DB{'gtdb'} =~ /gz$/){
    open DB,"gunzip -dc $DB{'gtdb'}|";
}else{
    open DB,"$DB{'gtdb'}";
}
my $head = <DB>;
chomp $head;
my @head = split /\t/,$head;
my %column = ();
for(my $i=0;$i<@head;$i++){
    $column{$head[$i]} = $i;
}

my @ass_list = ();
while(<DB>){
    chomp;
    my @l = split /\t/;
    if($l[$column{$search_id{$metaphlan_v}}] eq $taxid){
        push (@ass_list,$_);
    }
}

my $count = $#ass_list + 1;
print "Total genomes of $spe: $count\n";

my %select = ();
if($count < $threshold){
    for my $i (1..$count-1){
        $select{$i} = 1;
    }
}else{
    srand(999);
    for(1..$threshold){
        my $tmp;
        do {
            $tmp = int(rand($count));
        }while(exists $select{$tmp});
        $select{$tmp} = 1;
        #print "$tmp\n";
    }
}


foreach my $i(keys %select){
    my $ass_info = $ass_list[$i];
    my @ass_info = split /\t/,$ass_info;
    my ($id,$name) = ($ass_info[$column{'accession'}],$ass_info[$column{'ncbi_assembly_name'}]);
    my ($s1,$s2,$s3,$s4,$s5) = $id =~ /^.._((\w*)_(\d\d\d)(\d\d\d)(\d\d\d)\.\d)/;
    my $link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/$s2/$s3/$s4/$s5/$s1\_$name/$s1\_$name\_genomic.fna.gz";
    $link =~ s/ /_/g;
    print SH1 "wget -c -P $tmp_dir/dl/ $link\n";
}
print SH1 "gunzip $tmp_dir/dl/*gz\n";
print SH2 "for i in `ls $tmp_dir/dl/|grep fna|sed -e 's/_genomic.fna//'`;\ndo\n\tnucmer --mum -p $tmp_dir/snp/\$i $tmp_dir/$spe.marker.fa $tmp_dir/dl/\$i\\_genomic.fna;\n\tdelta-filter -1 $tmp_dir/snp/\$i.delta > $tmp_dir/snp/\$i.delta.2;\n\tshow-snps -CIrT $tmp_dir/snp/\$i.delta.2 > $tmp_dir/snp/\$i.snp;\ndone\n";


$DB{'fasttree'} = '/home/jiangyiqi/miniconda3/bin/fasttree';
open SH3,">$tmp_dir/3.tree.sh";
print SH3 "perl $RealBin/PT-06-1-narrow.pl $DB{'in'} $tmp_dir/snp/ $spe $tmp_dir/ $DB{'fasttree'}\n";
print SH3 "python $RealBin/nwk2mat.py $tmp_dir/tree.nwk $tmp_dir\n";




#------------------------------------#
#---------- sub-routines ------------#
#------------------------------------#

#---- test para and alert
sub para_alert{
	my $alert = 0;
    $tmp_dir = abs_path($tmp_dir);

	if($Help){
		$alert = 1;
	}
	elsif(!$tmp_dir || !-d $tmp_dir){
        `mkdir -p $tmp_dir`;    
	}else{
        print "Warning: working directory $tmp_dir exists.\n"
    }
    foreach my $s(@software){
        my $path = `which $s`;
        chomp $path;
        unless(-e $path){
            die "Check executeable of \"$s\".\n";
        }
    }
    foreach my $infile(keys %DB){
        if(!$DB{$infile}){
            $alert = 1;
        }
    }

	die "
	USAGE:
	perl $RealScript [<Options>] -WDR <whole_work_dir> -LS <strain_list> -PS <PStrain_dir> -S <similarity_cutoff>

	OPTIONS:
	-WDR  [s]  the whole work dir. <required>
	-S    [s]  the species name, e.g. s__Escherichia_coli. <required>
	-V    [s]  MetaPhlAn version M3/M4. <required>
	-I    [s]  PStrain output of strain sequence e.g. result/seq/s__Escherichia_coli_seq.txt. <required>
    -N    [s]  Number of randomly selected genomes. Default: 50
    -DBS  [s]  Species markers file of MetaPhlAn database. <required>
    -DBM  [s]  Marker genes fasta file from MetaPhlAn database. <required>
	-h         show this help.

	\n" if($alert == 1);

    foreach my $infile(keys %DB){
        unless(-e $DB{$infile}){
            die "Input not exists: $DB{$infile}.\n";
        }
        $DB{$infile} = abs_path($DB{$infile});
    }


	$_ = abs_path($_) for ($whole_dir);

}
