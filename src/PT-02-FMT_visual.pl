#!/usr/bin/perl  
#use strict;
use Cwd qw/abs_path/;
use Getopt::Long;
use FindBin qw/$RealBin $RealScript/;

my ($whole_dir,$experiment_list,$pstrain_dir,$similarity_cutoff);
my ($Help);

# embedded in the workflow
my $database_list = "$RealBin/config.list";

GetOptions(
	# required para
	"WDR=s"	=>	\$whole_dir,
	"LS=s"	=>	\$experiment_list,
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
my $vis_dir = "$whole_dir/FMT_visual";
`mkdir -p $vis_dir` if(!-d $vis_dir);


# link all run to lib
my %clu = ();
open (LIST,$experiment_list)||die"fail $experiment_list: $!\n";
while(<LIST>){
	chomp;
	my @sam = split /\t/;
	for(my $i=0;$i<@sam;$i++){
		my $type;
		if($i == 0){$type = 'donor'}elsif($i == 1){$type = 'receipt'}else{$type = 'postFMT'}
		my $sam = $sam[$i];
		open IN,"$pstrain_dir/$sam/result/strain_merged_RA_$similarity_cutoff.txt";
		<IN>;
		while(<IN>){
			chomp;
			my @l = split / /;
			$clu{$l[0]}{$type}{$l[3]} = 1;
		}
		close IN;
	}

}
close LIST;
my @order = ('donor','receipt','postFMT');
my %new = ();
foreach my $spe(sort keys %clu){
	my $n = 1;
	for(my $i =0;$i<@order;$i++){
		my %tmp = %{$clu{$spe}{$order[$i]}};
		foreach my $clu(keys %tmp){
			my $cluname = (sprintf "%02d", $n);
			unless(exists $new{$spe}{$clu}){
				$new{$spe}{$clu} = "clu$cluname";
				$n ++;
			}
		}
	}
}
open (LIST,$experiment_list)||die"fail $experiment_list: $!\n";
while(<LIST>){
	chomp;
	my @sam = split /\t/;
	my %all = ();
	for(my $i=0;$i<@sam;$i++){
		my $sam = $sam[$i];
		open IN,"$pstrain_dir/$sam/result/strain_merged_RA_$similarity_cutoff.txt";
		<IN>;
		while(<IN>){
			chomp;
			my @l = split / /;
			my $new = $new{$l[0]}{$l[3]};
			$all{"$l[0]_$new"}{$sam} += $l[5];
		}
		close IN;
	}
	my $h = join("\t",@sam);
	my $id = join("_",@sam);
	my $abd = "$vis_dir/$id.$similarity_cutoff.xls";
	open OU, ">$abd";
	print OU "tax\t$h\n";
	my %abd = ();
	my %rela = ();
	my %donor = ();
	my %recip = ();
	my %anno = ();
	foreach my $strain (sort keys %all){
		my $out = $strain;
		my ($spe,$clu) = split /_clu/,$strain;


		for(my $i=0;$i<@sam;$i++){
			my $sam = $sam[$i];
			if(exists $all{$strain}{$sam}){
				$out .= "\t$all{$strain}{$sam}";
				$abd{$spe}{$sam} .= "clu$clu:$all{$strain}{$sam};";
				$rela{$spe}{$sam} += $all{$strain}{$sam};

			}else{
				$out .= "\t0";
			}
		}
		print OU "$out\n";

		my @l = split /\t/,$out;
		my ($d,$r,$p) = ($l[1],$l[2],$l[3]);
		
		if($d > 0 && $p > 0 && $r == 0){
			$anno{$spe} .= "donor;";
		}elsif($d == 0 && $p > 0 && $r > 0){
			$anno{$spe} .= "receipt;";
		}elsif($p == 0 ){
			$anno{$spe} .= "not_colon;";
		}elsif($d >0 && $p > 0 && $r > 0){
			$anno{$spe} .= "same";
		}else{
			$anno{$spe} .= "uncertain;";
		}

	}
	close OU;

	my @order = ("donor","receipt","same","uncertain","not_colon");
	my %sort = ();
	my %out = ();
	foreach my $spe (keys %abd){
		my $out = $spe;
		my $not0 = 0;
		for(my $i=0;$i<3;$i++){
			my $sam = $sam[$i];
			if(exists $abd{$spe}{$sam}){
				$out .= "\t$abd{$spe}{$sam}";
				$not0 += 1 + 0.5*($#sam-$i);
			}else{
				$out .= "\t0";
			}
		}
		my @l = split /\t/,$out;
		my $anno = $anno{$spe};
		for(my $i=0;$i<@order;$i++){
			if($anno =~ /$order[$i]/){
				#my $type = $order[$i];
				$sort{$order[$i]}{$spe} = $not0;
				last;
			}
		}
		$out{$spe} = $out;
	}

	open OU,">$vis_dir/$id.$similarity_cutoff.rela.vis.txt";
	print OU "Species\tSample\tType\tSource\tAbd\n";
	my @type = @order;
	for($t=0;$t<@type;$t++){
		my %tmp = %{$sort{$type[$t]}};
		foreach my $spe(sort {$tmp{$b}<=>$tmp{$a}} keys %tmp){
			my @out = split /\t/,$out{$spe};
			for(my $i=0;$i<@sam;$i++){
				my $sam = $sam[$i];
				if(exists $abd{$spe}{$sam}){
					my @clu = split /;/,$abd{$spe}{$sam};
					foreach my $clu(@clu){
						my @clu_info = split /:/,$clu;
						my $rela = $clu_info[1]*100/$rela{$spe}{$sam};
						my $clu_name = $clu_info[0];
						print OU "$spe\t$sam\t$type[$t]\t$clu_name\t$rela\n";
					}
				
				}
			
			}
		}

	}
	close OU;
	my $order = join("','",@sam);
	open R,">$vis_dir/$id.$similarity_cutoff.rela.R";
	my $height = ($#sam+2)*1 + 1.5;
	$R = << "runR";
library(ggplot2,lib.loc="$DB{'R_lib'}")
data<-read.table("$vis_dir/$id.$similarity_cutoff.rela.vis.txt",sep = "\t",header = T)
data\$Species <- as.character(data\$Species)
data\$Species <- factor(data\$Species, levels=unique(data\$Species))
data\$Sample <- factor(data\$Sample, levels=c('$order'))
library(RColorBrewer)
colourCount <- length(unique(data\$Source))

pdf('$vis_dir/$id.$similarity_cutoff.rela.pdf',width=12,height=$height)

g<-ggplot(data=data, aes(x=Species,fill=Source,y=Abd)) + geom_bar(position="stack",stat='identity',color="black") +facet_grid(Sample~.) + labs(title="Strain proportion change in FMT of $sam[1]", x ="Species", y = "proportion") +theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.9,size = 7,face = "italic"))
g +  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "RdYlBu"))(colourCount))

dev.off()

runR
print R $R;


system("$DB{'Rscript'}  $vis_dir/$id.$similarity_cutoff.rela.R");

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
	-h         show this help.

	\n" if($alert);

	$_ = abs_path($_) for ($whole_dir,$experiment_list,$pstrain_dir);

}
