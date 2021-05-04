#!/usr/bin/perl  
#use strict;
use Cwd qw/abs_path/;
use Getopt::Long;
use FindBin qw/$RealBin $RealScript/;

my ($whole_dir,$group_info,$pstrain_dir,$similarity_cutoff,$name);
my ($Help);

# embedded in the workflow
my $database_list = "$RealBin/config.list";
$name = 'output';

GetOptions(
	# required para
	"WDR=s"	=>	\$whole_dir,
	"LS=s"	=>	\$group_info,
	"PS=s"  =>	\$pstrain_dir,
	"S=s"	=>	\$similarity_cutoff,
	"PRE=s" =>	\$name,
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
my $test_dir = "$whole_dir/wilcox_testing";
`mkdir -p $test_dir` if(!-d $test_dir);


# link all run to lib
my %all = ();
my %spe = ();
my @sam = ();
open (LIST,$group_info)||die"fail $group_info: $!\n";
`cp $group_info $test_dir/$name.$similarity_cutoff.group_info`;
<LIST>;
while(<LIST>){
	chomp;
	my ($sam,$g) = split /\t/;
	push (@sam,$sam);
	open IN,"$pstrain_dir/$sam/result/strain_merged_RA_$similarity_cutoff.txt";
	<IN>;
	while(<IN>){
		chomp;
		my @l = split / /;
		$all{"$l[0]_$l[3]"}{$sam} += $l[5];
		$spe{"$l[0]"}{$sam} += $l[5];
	}
}
close LIST;


open OU,">$test_dir/$name.strain.$similarity_cutoff.tsv";
open OT,">$test_dir/$name.species.$similarity_cutoff.tsv";

@sam = sort @sam;
my $h = join("\t",@sam);
print OU "ID\t$h\n";
print OT "ID\t$h\n";

my $cutoff = 0.2;
foreach my $strain (sort keys %all){
	my $out = $strain;
	my $not0 = 0;
	for($i=0;$i<@sam;$i++){
		my $sam = $sam[$i];
		if(exists $all{$strain}{$sam}){
			$out .= "\t$all{$strain}{$sam}";
			$not0 ++;
		}else{
			$out .= "\t0";
		}
	}
	if($not0 > ($#sam+1)*$cutoff){
		print OU "$out\n";
	}
}

foreach my $spe (sort keys %spe){
	my $out = $spe;
	my $not0 = 0;
	for($i=0;$i<@sam;$i++){
		my $sam = $sam[$i];
		if(exists $spe{$spe}{$sam}){
			$out .= "\t$spe{$spe}{$sam}";
			$not0 ++;
		}else{
			$out .= "\t0";
		}
	}
	if($not0 > ($#sam+1)*$cutoff){
		print OT "$out\n";
	}
}




open R ,">$test_dir/$name.$similarity_cutoff.R";
$R = << "runR";
infile <- file("$test_dir/$name.strain.$similarity_cutoff.tsv","r")
outfile <- file("$test_dir/$name.strain.$similarity_cutoff.test.tsv","w")
line <- readLines(infile, n=1)
samp <- as.vector(unlist(strsplit(line, split="\t")))
phen <- read.table("$test_dir/$name.$similarity_cutoff.group_info",head=T)
phen_samp <- cbind("samp"=as.vector(phen[,1]), "state"=as.vector(phen[,2]))
groupname <- as.matrix(sort(unique(phen\$groupname)))
title <- "ID"
for (i in 1:nrow(groupname))
{
	mean <- paste("mean(",groupname[i,1],")",sep="")
	sd <- paste("sd(",groupname[i,1],")",sep="")
	occ <- paste("occ-rate(",groupname[i,1],")",sep="")
	title <- paste(title,mean,sd,occ,sep="\t")
}
outmatrixname <- paste(title,"enriched","pvalue",sep="\t")
writeLines(outmatrixname, con=outfile, sep="\n")


while(length(line <- readLines(infile, n=1))) {
	line <- as.vector(unlist(strsplit(line, split="\t")))
	tmp <- cbind(samp=samp[-1], abund=line[-1])
	dat <- merge(phen_samp, tmp, by="samp",sort =F)
	dat\$abund <- as.numeric(as.vector((dat\$abund)))
	dat <- dat[complete.cases(dat),]
	mean.sd <- as.matrix(aggregate(dat\$abund,by=list(dat\$state),FUN=function(x)c(mean=sprintf("%0.9f",mean(x)),sd=sprintf("%0.9f",sd(x)))))
	testresult <- wilcox.test(dat\$abund ~ dat\$state,alternative ="two.sided")
	test.pvalue <- testresult\$p.value
	Rank <- rank(dat\$abund)
	Rank_mean <- as.matrix(aggregate(Rank,by=list(dat\$state),FUN=function(x)c(mean=mean(x))))
	Occ <- as.matrix(aggregate(dat\$abund,by=list(dat\$state),FUN=function(x)sum(x != 0)/length(x)))
	enriched <- as.character(Rank_mean[which.max(Rank_mean[,2]),1])
	output <- paste(line[1],mean.sd[1,2],mean.sd[1,3],format(as.numeric(Occ[1,2]),digits=3),mean.sd[2,2],mean.sd[2,3],format(as.numeric(Occ[2,2]),digits=3),enriched,test.pvalue,sep="\t")
	writeLines(output, con=outfile, sep="\n")
}
close(infile)
close(outfile)

data <- read.table("$test_dir/$name.strain.$similarity_cutoff.test.tsv",head=T,check.names=F,sep="\t")
p_adjusted <- p.adjust(data\$pvalue, method="BH")
data\$qvalue <- p_adjusted
write.table(format(data,scientific=FALSE,digits=9),file="$test_dir/$name.strain.$similarity_cutoff.test.tsv",row.names=F, col.names=T, quote=F, sep="\t")
runR
print R $R;

my $R2 = $R;
$R2 =~ s/$name\.strain\.$similarity_cutoff/$name.species.$similarity_cutoff/g;
print R $R2;

system("$DB{'Rscript'}  $test_dir/$name.$similarity_cutoff.R");

my $p_cutoff = 0.05;
my $test1 = "$test_dir/$name.strain.$similarity_cutoff.test.tsv"; 
my $test2 = "$test_dir/$name.species.$similarity_cutoff.test.tsv";
my $abd = "$test_dir/$name.strain.$similarity_cutoff.tsv";
%pval = ();
my @file = ($test1,$test2);
foreach my $test(@file){
	open IN,$test;
	<IN>;
	while(<IN>){
		chomp;
		my @l = split /\t/;
		if($l[-1] < $p_cutoff){
			$pval{$l[0]} = 'S';
		}else{
			$pval{$l[0]} = 'NS';
		}
	}
}

open IN,$abd;
open OU,">$test_dir/$name.$similarity_cutoff.heatmap.tsv";
my $h = <IN>;
@h = split /\t/,$h;
$h[0] = "$h[0]\tStrain\tSpecies";
$h = join("\t",@h);
my $n = 0;
print OU $h;
while(<IN>){
	chomp;
	my @l = split /\t/;
	my $spe = (split /_clu/,$l[0])[0];
	$l[0] = "$l[0]\t$pval{$l[0]}\t$pval{$spe}";
	$l = join("\t",@l);
	print OU "$l\n";
	$n ++;

}

my $height = 3 + $n*0.2;
open R,">$test_dir/$name.$similarity_cutoff.heatmap.R";
$R = << "runR";
library(pheatmap,lib.loc="$DB{'R_lib'}")
data <- read.table("$test_dir/$name.$similarity_cutoff.heatmap.tsv",header=1,row.names = 1,sep="\t")
input <- data[,-c(1,2)]
anno <- data[,c(1,2)]
anno2 <- read.table("$test_dir/$name.$similarity_cutoff.group_info",header=1,row.names = 1,sep="\t")
colnames(anno2) <- c("Group")
input<-log10(input+0.00001)

pdf('$test_dir/$name.$similarity_cutoff.heatmap.pdf',width=10,height=$height)
pheatmap(input,annotation_row = anno,cutree_cols = 2,border=T,annotation_col = anno2,main = "-log10(relative abundance)")
dev.off()

runR
print R $R;

system("$DB{'Rscript'}  $test_dir/$name.$similarity_cutoff.heatmap.R");


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
	elsif(!$group_info || !-e $group_info){
		$alert = 1;
	}
	elsif(!$pstrain_dir || !-d $pstrain_dir){
		$alert = 1;
	}


	die "
	USAGE:
	perl $RealScript [<Options>] -WDR <whole_work_dir> -LS <sample_group_list> -PS <PStrain_dir> -S <similarity_cutoff> 

	OPTIONS:
	-WDR  [s]  the whole work dir. <required>
	-LS   [s]  the sample list with group info. One sample a line, sample and the group separated by tab. Header line is needed <required>
	-PS   [s]  PStrain result direction. <required>
	-S    [s]  similarity cutoff. <required>
	-PRE  [s]  prefix of output. [default=output]
	-h         show this help.

	\n" if($alert);

	$_ = abs_path($_) for ($whole_dir,$group_info,$pstrain_dir);

}
