#!/usr/bin/perl -w
use strict;

#my ($clu,$snp,$spe,$out,$fasttree) = ('demo_PStrain_result//merge_0.9/seq/Ruminococcus_lactaris_clu.txt','find_sim/Ruminococcus_lactaris_clu-1/snp/','Ruminococcus_lactaris','find_sim/Ruminococcus_lactaris_clu-1/','/mnt/disk2_workspace/jiangyiqi/bin/FastTree');
my ($clu,$snp,$spe,$out,$fasttree) = @ARGV;
open IN,$clu;
my $h = <IN>;
chomp $h;
$h =~ s/\//_/g;
my @h = split /\t/,$h;
my %in = ();
my %all = ();
while(<IN>){
	chomp;
	my @l = split /\t/;
	my ($ref,$alt) = ($l[2],$l[3]);
	for(my $i=4;$i<@l;$i++){
		$in{"$l[0]_$l[1]"} = $ref;
		if($l[$i] eq 1){
			$all{"$spe\_$h[$i]"}{"$l[0]_$l[1]"} = $alt;
		}else{
			$all{"$spe\_$h[$i]"}{"$l[0]_$l[1]"} = $ref;
		}
	}
}
my @file = glob("$snp/*snp");
foreach my $file(@file){
	open IN,$file;
	$file =~ /$snp\/(.*)\.snp/;
	my $id = $1;
	<IN>;<IN>;<IN>;<IN>;
	while(<IN>){
		chomp;
		my @l = split /\t/;
		#if(!exists $in{"$l[8]_$l[0]"}){next}
		$in{"$l[8]_$l[0]"} = 1;
		if(length($l[2])>1){next}
		$all{$id}{"$l[8]_$l[0]"} = $l[2];
	}
}

open OU,">$out/tree.fa";
foreach my $sam (keys %all){
	my $out = ">$sam\n";
	foreach my $pos(sort keys %in){
		if(exists $all{$sam}{$pos}){
			$out .= $all{$sam}{$pos};
		}else{
			$out .= $in{$pos};
		}
	}
	print OU "$out\n";
}
system("$fasttree -nt $out/tree.fa > $out/tree.nwk");
