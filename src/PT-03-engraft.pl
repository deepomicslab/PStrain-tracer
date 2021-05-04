#!/usr/bin/perl -w 
use strict;
use Cwd qw/abs_path/;
use Getopt::Long;
use FindBin qw/$RealBin $RealScript/;

my ($whole_dir,$experiment_list,$pstrain_dir,$similarity_cutoff,$mode);
my ($Help);

# embedded in the workflow
my $database_list = "$RealBin/config.list";
$mode = 'propotion';
GetOptions(
	# required para
	"WDR=s"	=>	\$whole_dir,
	"LS=s"	=>	\$experiment_list,
	"PS=s"  =>	\$pstrain_dir,
	"S=s"	=>	\$similarity_cutoff,
	"M=s"	=>	\$mode,
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
my $vis_dir = "$whole_dir/FMT_visual";
my $en_dir = "$whole_dir/engraft";
die "no $vis_dir found! should run PT-02 first.\n" if(!-d $vis_dir);
`mkdir -p $en_dir` if(!-d $en_dir);


# link all run to lib
my %clu = ();
open (LIST,$experiment_list)||die"fail $experiment_list: $!\n";
open OU,">$en_dir/engraft.$similarity_cutoff.$mode.tsv";
print OU "ABD\tGROUP\n";
while(<LIST>){
	chomp;
	my @sam = split /\t/;
	my $fid = join("_",@sam);
	my $file = "$vis_dir/$fid.$similarity_cutoff.xls";
	open IN,$file;
	<IN>;
	my (%in,%n,%total) = ();
	while(<IN>){
		chomp;
		my @l = split /\t/;
		my @info = split /_clu/,$l[0];
		if($l[2] > 0 && $l[3] > 0){$in{$info[0]} = 1;}
		if($l[2] > 0){$n{$info[0]} ++ }
		$total{$info[0]} += $l[2];
	}

	close IN;
	open IN,$file;
	<IN>;
	while(<IN>){
		chomp;
		my @l = split /\t/;
		my $per = '';
		my @info = split /_clu/,$l[0];
		if(!exists $in{$info[0]}||$n{$info[0]} eq 1){next}
		if($l[2] > 0){
			if($mode eq 'abd'){
				$per = $l[2];
			}else{
				$per = $l[2]/$total{$info[0]};
			}
			if($l[3] > 0){
				print OU "$per\tengraft\n";
			}else{
				print OU "$per\tno_engraft\n";
			}

		}
		
	}
}
close LIST;

open R,">$en_dir/engraft.$similarity_cutoff.$mode.R";
my $R = << "runR";
library(ggplot2,lib.loc="$DB{'R_lib'}")
data<-read.table("$en_dir/engraft.$similarity_cutoff.$mode.tsv",sep = "\t",header = T)
x<-subset(data,GROUP == "engraft")\$ABD
y<-subset(data,GROUP == "no_engraft")\$ABD
p<-wilcox.test(x,y)\$p.value
lab<-sprintf("p = %.2f",p)
pdf('$en_dir/engraft.$similarity_cutoff.$mode.pdf',width=4,height=5)

ggplot(data, aes(x=GROUP, y=ABD,fill=GROUP)) + geom_boxplot()+ theme_bw()+ annotate("text", x=1.5, y=max(data\$ABD)*1.1, label= lab)+xlab("Patient")+ylab("$mode")+ theme(legend.position = "none")

dev.off()

runR
print R $R;


system("$DB{'Rscript'} $en_dir/engraft.$similarity_cutoff.$mode.R");


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
	elsif(!$experiment_list || !-e $experiment_list){
		$alert = 1;
	}
	elsif(!$pstrain_dir || !-d $pstrain_dir){
		$alert = 1;
	}

	die "
	USAGE:
	perl $RealScript [<Options>] -WDR <whole_work_dir> -LS <experiment_list> -PS <PStrain_dir> -S <similarity_cutoff>

	OPTIONS:
	-WDR  [s]  the whole work dir. <required>
	-LS   [s]  the sample-list of each FMT experiment in one line, separated by tab. <required>
	-PS   [s]  PStrain result direction. <required>
	-S    [s]  similarity cutoff. <required>
	-M    [s]  the acceptable input: 'abd' or 'propotion'. The 'abd' means the origin abundance,
       	           the 'propotion' means the normalized abundance. [default:propotion]
	-h         show this help.

	\n" if($alert);

	$_ = abs_path($_) for ($whole_dir,$experiment_list,$pstrain_dir);

}
