use strict;
use warnings;
use File::Basename qw(basename dirname);



my $sequence_file=$ARGV[0];
my $David_cluster_file=$ARGV[1];
my $gene_info_file=$ARGV[2];
my $calling_gene_file=$ARGV[3];
my $half_length=$ARGV[4];
my $summit=$ARGV[5];


$half_length=100 if(!defined $half_length);
$summit=3 if(!defined $summit);
$summit--;

if(!-d "$David_cluster_file.meme"){
	system("mkdir $David_cluster_file.meme");
}

print "$sequence_file\n$David_cluster_file\n$gene_info_file\n";

my %ID2gene=();
open(IN,$gene_info_file)||die;
<IN>;
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	$ID2gene{$a[3]}{$a[8]}=0;
}
close IN;

my %David_cluster=();
open(IN,$David_cluster_file)||die;
my $na="";
my $es="";
while(my $t=<IN>){
	next if($t=~/^$/);
	next if($t=~/^Category/);
	if($t=~/^Annotation Cluster (\d+).*Enrichment Score: (.*?)$/){
		$na=$1;
		$es=round($2,8);
	}
	else{
		chomp $t;
		my @a=split(/\t/,$t);
		my @genes=split(/, /,$a[5]);
		foreach my $id (@genes){
			foreach my $g (keys %{$ID2gene{$id}}){
				$David_cluster{("C$na"."_$es")}{$g}=0;
			}
		}
	}
}
close IN;

my %sequence=();
my $loc="";
open(IN,$sequence_file)||die;
while(my $t=<IN>){
	chomp $t;
	if($t=~/^>(.*)$/){
		$loc=$1;
	}
	else{
		$sequence{$loc}.=$t;
	}
}
close IN;

open(IN,$calling_gene_file)||die;
my %calling_gene=();
my $str=<IN>;
chomp $str;
my @header=split(/\t/,$str);

my $gene_body="";
my $first_choice_gene="";
my $first_choice_distance="";
my $second_choice_gene="";
my $second_choice_distance="";
my $ABCD="";
my $modified_upstream_target="";
my $modified_upstream_target_dist="";
my $modified_dnstream_target="";
my $modified_dnstream_target_dist="";
my $i=0;
foreach my $t (@header){
	$gene_body=$i if($t =~/^gene_body/);
	$first_choice_gene=$i if($t eq "first_choice_gene");
	$first_choice_distance=$i if($t eq "first_choice_distance");
	$second_choice_gene=$i if($t eq "second_choice_gene");
	$second_choice_distance=$i if($t eq "second_choice_distance");
	$ABCD=$i if($t eq "ABCD");
	$modified_upstream_target=$i if($t eq "modified_upstream_target");
	$modified_upstream_target_dist=$i if($t eq "modified_upstream_target_dist");
	$modified_dnstream_target=$i if($t eq "modified_dnstream_target");
	$modified_dnstream_target_dist=$i if($t eq "modified_dnstream_target_dist");
	$i++;
}

while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	if($a[$ABCD] eq "A" or $a[$ABCD] eq "D" or $a[$ABCD] eq "C"){
		my $g=$a[$first_choice_gene];
		$g=~s/,.*//;
		my $loc="$a[0]:".($a[$summit]-100)."-".($a[$summit]+100);
		$calling_gene{$g}{$loc}=0;
	}
	else{
		if($a[$modified_upstream_target_dist] ne "-" and $a[$modified_upstream_target_dist]<2000){
			my $g=$a[$modified_upstream_target];
			$g=~s/,.*//;
			my $loc="$a[0]:".($a[$summit]-100)."-".($a[$summit]+100);
			$calling_gene{$g}{$loc}=0;
		}
		if($a[$modified_dnstream_target_dist] ne "-" and $a[$modified_dnstream_target_dist]<2000){
			my $g=$a[$modified_dnstream_target];
			$g=~s/,.*//;
			my $loc="$a[0]:".($a[$summit]-100)."-".($a[$summit]+100);
			$calling_gene{$g}{$loc}=0;
		}
		
	}
	
	if($a[$ABCD] eq "C"){
		if($a[$second_choice_distance] ne "-" and $a[$second_choice_distance]<2000){
			my $g=$a[$second_choice_gene];
			$g=~s/,.*//;
			my $loc="$a[0]:".($a[$summit]-100)."-".($a[$summit]+100);
			$calling_gene{$g}{$loc}=0;
		}
	}

}
close IN;

foreach my $c (sort keys %David_cluster){
	open(OUT,">$David_cluster_file.meme/$c.fasta")||die;
	my $line_count=0;
	foreach my $g (sort keys %{$David_cluster{$c}}){
		foreach my $loc (sort keys %{$calling_gene{$g}}){
			print OUT ">$loc $g\n$sequence{$loc}\n";
			$line_count++;
		}
	}
	close OUT;
	my $command="";
	$command.="meme-chip $David_cluster_file.meme/$c.fasta ";
	$command.="-oc $David_cluster_file.meme/$c.memechip ";
	$command.="-ccut 0 -meme-mod zoops -meme-minw 6 -meme-maxw 15 -meme-nmotifs 10 -meme-maxsize ".($line_count*250)." ";
	$command.="-db /Users/zqmiao/Dropbox\\ \\(Personal\\)/mytools/motif_databases/JASPAR/JASPAR2018_CORE_fungi_non-redundant.meme";
	
	print $command."\n";
	system($command);
	
}





sub round{
	my ($s,$c)=@_;#$c, how many decimals want to keep
	my $n1=10**$c;
	#print "$c\t$n1\n";
	return $s=int(($s*$n1)+0.5)/$n1;
}







