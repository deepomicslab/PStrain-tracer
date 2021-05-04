#!/usr/bin/perl  
#use strict;
use Cwd qw/abs_path/;
use Getopt::Long;
use FindBin qw/$RealBin $RealScript/;

my ($whole_dir,$pstrain_dir,$similarity_cutoff);
my ($Help);

# embedded in the workflow
my $database_list = "$RealBin/config.list";

GetOptions(
	# required para
	"WDR=s"	=>	\$whole_dir,
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
my $strN_dir = "$whole_dir/strainN";
`mkdir -p $strN_dir` if(!-d $strN_dir);


# link all run to lib

open IN,"$pstrain_dir/merge_$similarity_cutoff/strain_number.txt" || die "need run PT-01 first: $!\n";
my $h = <IN>;
chomp $h;
my @h = split /\t/,$h;
my %all = ();
while(<IN>){
	chomp;
	my @l = split /\t/;
	for(my $i=1;$i<@l;$i++){
		if($l[$i] ne 0){
			$all{$l[0]}{$h[$i]} = $l[$i];
		}
	}
}

open OU,">$strN_dir/strainN_heterSNV_stat.$similarity_cutoff.tsv";
print OU "Sample\tSpecies\tStrain_N\tHeterzygous_SNV_N\n";
foreach my $sam(keys %all){
	my @spe = keys %{$all{$sam}};
	open IN,"gunzip -dc $pstrain_dir/$sam/map/mapped.vcf.gz|";
	my %n = ();
	while(<IN>){
		if(/^#/){next}
		chomp;
		foreach my $spe(@spe){
			if($_ =~ /$spe/){
				chomp;
				my @l = split /\t/;
				my $gt = (split /:/,$l[9])[0];
				if($gt eq '0/1'){$n{$spe} ++}
			}
		}
	}
	foreach my $spe(@spe){
		my $n = 0;
		if(exists $n{$spe}){
			$n = $n{$spe};
		}
		print OU "$sam\t$spe\t$all{$sam}{$spe}\t$n\n";
	}
}




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
	elsif(!$pstrain_dir || !-d $pstrain_dir){
		$alert = 1;
	}

	die "
	USAGE:
	perl $RealScript [<Options>] -WDR <whole_work_dir>  -PS <PStrain_dir> -S <similarity_cutoff>

	OPTIONS:
	-WDR  [s]  the whole work dir. <required>
	-PS   [s]  PStrain result direction. <required>
	-S    [s]  similarity cutoff. <required>
	-h         show this help.

	\n" if($alert);

	$_ = abs_path($_) for ($whole_dir,$pstrain_dir);

}
