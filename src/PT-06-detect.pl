#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use Getopt::Long;
use FindBin qw/$RealBin $RealScript/;

my ($whole_dir,$stain_list,$pstrain_dir,$similarity_cutoff);
my ($Help);

# embedded in the workflow
my $database_list = "$RealBin/config.list";

GetOptions(
	# required para
	"WDR=s"	=>	\$whole_dir,
	"LS=s"	=>	\$stain_list,
	"PS=s"  =>	\$pstrain_dir,
	"S=s"	=>	\$similarity_cutoff,
	# help
	"h"		=>	\$Help
);

# check para
&para_alert;

#----------- read database list -----------
my %DB;
open(DBL,$database_list)||die"fail $database_list: $!\n";
while(<DBL>){
	next if(/^\s*\#/ || /^\s+$/);
	s/\s//g;
	my ($key,$value) = (/^([^\=]+)\=([^\=]+)$/);
	$DB{$key} = $value;
}
close DBL;

#----------- dir --------
my $detect_dir = "$whole_dir/find_strain";
`mkdir -p $detect_dir` if(!-d $detect_dir);


# link all run to lib
open (LIST,$stain_list)||die"fail $stain_list: $!\n";
while(my $str = <LIST>){
	chomp $str;
	my $spe = (split /_clu-/,$str)[0];
	my $tmp_dir = "$detect_dir/$str";
	`mkdir -p $tmp_dir`;
	`zcat $DB{'species_marker'}|grep $spe > $tmp_dir/$spe.list`;
	open IN,"$tmp_dir/$spe.list";
	my %in = ();
	while(<IN>){
		chomp;
		my @l = split /\s/;
		$in{$l[0]} = 1;
	}
	open OU,">$tmp_dir/$spe.marker.fa";
	$/ = '>';
	open IN,"gunzip -dc $DB{'marker_gene'}|";
	<IN>;
	while(my $seq=<IN>){
		$seq =~ s/>$//;
		my $id = (split /\n/,$seq)[0];
		if(exists $in{$id}){
			print OU ">$seq";
		}
	}
	$/ = "\n";
	open SH1,">$tmp_dir/1.dl.sh";
	`mkdir -p $tmp_dir/dl/`;
	
	open SH2,">$tmp_dir/2.snp.sh";
	`mkdir -p $tmp_dir/snp/`;

	my @spe_info = split /_/,$spe;
	my $ass_list = `cat $DB{'ass_sum'}|grep $spe_info[0]|grep $spe_info[1]`;
	my @ass_list = split /\n/,$ass_list;
	
	foreach my $ass_info(@ass_list){
		chomp $ass_info;
		my @ass_info = split /\t/,$ass_info;
		my $link = "$ass_info[19]/$ass_info[0]_$ass_info[15]_genomic.fna.gz";
		$link =~ s/ftp/http/;
		print SH1 "wget -c -P $tmp_dir/dl/ $link\n";
	}
	print SH1 "gunzip $tmp_dir/dl/*gz\n";
	print SH2 "for i in `ls $tmp_dir/dl/|grep fna|sed -e 's/_genomic.fna//'`;\ndo\n\t$DB{'MUMmer'}/nucmer --mum -p $tmp_dir/snp/\$i $tmp_dir/$spe.marker.fa $tmp_dir/dl/\$i\\_genomic.fna;\n\t$DB{'MUMmer'}/delta-filter -1 $tmp_dir/snp/\$i.delta > $tmp_dir/snp/\$i.delta.2;\n\t$DB{'MUMmer'}/show-snps -CIrT $tmp_dir/snp/\$i.delta.2 > $tmp_dir/snp/\$i.snp;\ndone\n";
	
	open SH3,">$tmp_dir/3.tree.sh";
	print SH3 "perl $RealBin/PT-06-1-narrow.pl $pstrain_dir/merge_$similarity_cutoff/seq/$spe\_clu.txt $tmp_dir/snp/ $spe $tmp_dir/ $DB{'fasttree'}\n";

}
close LIST;



#------------------------------------#
#---------- sub-routines ------------#
#------------------------------------#

#---- test para and alert
sub para_alert{
	my $alert;

	if($Help){
		$alert = 1;
	}
	elsif(!$whole_dir || !-d $whole_dir){
		$alert = 1;
	}
	elsif(!$stain_list || !-e $stain_list){
		$alert = 1;
	}
	elsif(!$pstrain_dir || !-d $pstrain_dir){
		$alert = 1;
	}


	die "
	USAGE:
	perl $RealScript [<Options>] -WDR <whole_work_dir> -LS <strain_list> -PS <PStrain_dir> -S <similarity_cutoff>

	OPTIONS:
	-WDR  [s]  the whole work dir. <required>
	-LS   [s]  the strain list. <required>
	-PS   [s]  PStrain result direction. <required>
	-S    [s]  similarity cutoff. <required>
	-h         show this help.

	\n" if($alert);

	$_ = abs_path($_) for ($whole_dir,$stain_list,$pstrain_dir);

}
