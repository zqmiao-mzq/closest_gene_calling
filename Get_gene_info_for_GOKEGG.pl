use strict;
use warnings;
#use LWP::Simple;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use URI::Escape;
use Pod::Usage;

###############  Start Time  ########################################################
my $Time_Start="";
my $TS=time();
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";

##################  Main     ######################################################

#my $input1="./L1vsE1_up";
#my $input2="./L1vsE1_dn";


#perl ../../mytools/perl/Get_gene_info_for_GOKEGG.pl re_group_hsf1_targets.list.txt AN

my $input_file=$ARGV[0];

#my $input1=$ARGV[0];
#my $input2=$ARGV[1];
my $specis=$ARGV[1]; #AN,AF,CA,CG,Mouse,Human,ZebraFish
#my $categories = $ARGV[2];
#my $categories = "GOTERM_BP_4,GOTERM_MF_4,GOTERM_BP_FAT,GOTERM_MF_FAT,KEGG";
my $categories = "GOTERM_BP_4,GOTERM_MF_4,GOTERM_BP_5,GOTERM_MF_5,GOTERM_BP_FAT,GOTERM_MF_FAT,GOTERM_BP_DIRECT,GOTERM_MF_DIRECT,KEGG";
#my $categories = "GOTERM_BP_FAT";


my $gene_list="/Users/zqmiao/Dropbox\ \(Personal\)/Sties/chris_lab/KEGG_GO_dir/Gene_INFO/All_$specis.gene_info";
my $chartReport_file="/Users/zqmiao/Dropbox\\ \\(Personal\\)/mytools/perl/chartReport_for_GO.pl";
my $KEGG_List_file="/Users/zqmiao/Dropbox\ \(Personal\)/Sties/chris_lab/KEGG_GO_dir/KEGG_GENE/All_$specis.KEGG";
my $KEGG_GENE_file="/Users/zqmiao/Dropbox\ \(Personal\)/Sties/chris_lab/KEGG_GO_dir/KEGG_GENE/$specis.KEGG_GENE";


my %NCBI_ID=();
my %NCBI_Synonyms_ID=();
my %Info=();
my %NCBI_ID2KEGG_ID=();
open(IN,"$gene_list")||die "can't read $gene_list\n";
<IN>;
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	$NCBI_ID2KEGG_ID{$a[1]}=$a[5];
	$NCBI_ID{$a[1]}=$a[1];
	$Info{$a[1]}=$t;
	$NCBI_ID{uc($a[2])}=$a[1];
	$NCBI_ID{$a[3]}=$a[1];
	$NCBI_ID{$a[6]}=$a[1];
	$NCBI_ID{$a[7]}=$a[1];
	
}
close IN;

open(IN,"$gene_list")||die;
<IN>;
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	if($a[4] ne "-"){
		my @aa=split(/\|/,$a[4]);
		foreach my $tt (@aa){
			$NCBI_ID{$tt}=$a[1] if(!exists $NCBI_ID{$tt});
		}
	}
}
delete $NCBI_ID{"-"};
close IN;

# open(IN,"$input_file")||die "can't read $input_file\n";
# #open(OUT,">$input_file.test")||die "can't read $input_file\n";
# my $str=<IN>;
# chomp $str;
# my @group_name=split(/\t/,$str);
# my %input_mat=();
# while(my $t=<IN>){
# 	chomp $t;
# 	my @a=split(/\t/,$t);
# 	for(my $i=0;$i<@a;$i++){
# 		push(@{$input_mat{$group_name[$i]}},$a[$i]) if($a[$i] ne "");
# 		#print OUT "$a[$i]\t";
# 	}
# 	#print OUT "\n";
# }
# close IN;
# #close OUT;


open(IN,"$input_file")||die "can't read $input_file\n";
#open(OUT,">$input_file.test")||die "can't read $input_file\n";

my %input_mat=();
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	#push(@{$input_mat{$a[1]}},$a[0]);

	for(my $i=1; $i<@a;$i++){
		push(@{$input_mat{$a[0]}},$a[$i]);
	}
}
close IN;
#close OUT;
my @group_name=sort keys %input_mat;

print "@group_name\n";



my $path=dirname($input_file);

open(OUT,">$input_file.gene_info")||die "$input_file.gene_info";
print OUT "Group\tInput1\tTax_ID\tNCBI_ID\tSymbol\tLocusTag\tSynonyms\tKEGG_ID\tLab_ID\tCommonName\tdescription";

my %entry_KEGG=();
my %entry_GO=();


foreach my $group (@group_name){
	foreach my $g (@{$input_mat{$group}}){
		if(exists $NCBI_ID{$g}){
			push(@{$entry_GO{$group}},$NCBI_ID{$g});
			push(@{$entry_KEGG{$group}},$NCBI_ID{$g}) if(exists $NCBI_ID2KEGG_ID{$NCBI_ID{$g}});
			print OUT "$group\t$g\t$Info{$NCBI_ID{$g}}\n";
		}
		else{
			print OUT "$group\t$g\t-\t-\t-\t-\t-\t-\t-\t-\t-\t\n";
		}
	}
}

close OUT;


####  RUN GO  

open(OUT,">$input_file.GO_David")||die;
#open(OUT1,">$input_file.run_go.command")||die;

my %input_mat_gene_count=();

foreach my $group (@group_name){
	print "######$group\n";
	my $n=@{$entry_GO{$group}};
	$input_mat_gene_count{$group}=$n;
	my $input_gene_list=join(",",@{$entry_GO{$group}});
	system("perl $chartReport_file $input_gene_list $categories $input_file.GO_David $group");
	#print OUT1 "perl $chartReport_file $input_gene_list $categories $input_file.GO_David $group\n";
}

close OUT;
#close OUT1;
#####  RUN KEGG  


open(IN,$KEGG_List_file)||die "can't read $KEGG_List_file\n";
my %kegg_name=();
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	$kegg_name{$a[0]}=$a[1];
	
}
close IN;
open(IN,$KEGG_GENE_file)||die "can't read $KEGG_GENE_file\n";
my %kegg_gene=();
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	$kegg_gene{$a[1]}{$a[2]}=0;
}
close IN;


my %kegg_gene_mat=();
foreach my $group (@group_name){
	if(exists $entry_KEGG{$group}){
		foreach my $k (sort keys %kegg_name){
			foreach my $g (@{$entry_KEGG{$group}}){
				if(exists $kegg_gene{$k}{$NCBI_ID2KEGG_ID{$g}}){
					push(@{$kegg_gene_mat{$k}{$group}},$g);
				}
			}
		}
	}
}
open(OUT,">$input_file.kegg_mat")||die;
open(OUT1,">$input_file.kegg_gene_mat")||die;

print OUT "KEGG\tName";
print OUT1 "KEGG\tName";
foreach my $group (@group_name){
	print OUT "\t$group";
	print OUT1 "\t$group";
}
print OUT "\tall_pathway_gene_count\n";

print OUT1 "\n";

foreach my $k (sort keys %kegg_gene_mat){
	print OUT "$k\t$kegg_name{$k}";
	print OUT1 "$k\t$kegg_name{$k}";
	my @count_arr=();
	my $n=keys %{$kegg_gene{$k}};
	foreach my $group (@group_name){
		if(exists $kegg_gene_mat{$k}{$group}){
			my $c=@{$kegg_gene_mat{$k}{$group}};
			push(@count_arr,$c);
			print OUT1 "\t".join(",",@{$kegg_gene_mat{$k}{$group}});
			print OUT "\t$c";
		}
		else{
			push(@count_arr,0);
			print OUT1 "\t-";
			print OUT "\t0";
		}
	}
	print OUT "\t$n\n";
	
	print OUT1 "\n";
	
}

print OUT "all_count_by_group\tall_count_by_group";
foreach my $group (@group_name){
	print OUT "\t$input_mat_gene_count{$group}";
}
my $all_gene_count=keys %NCBI_ID;
print OUT "\t$all_gene_count\n";

close OUT;
close OUT1;

### summary GO mat


open(IN,"$input_file.GO_David")||die;
my %GO_David_mat=();
while(my $t=<IN>){
	next if($t=~/^Group/);
	chomp $t;
	my @a=split(/\t/,$t);
	while($a[2]=~/'/){
		$a[2]=~s/'//;
	}
	$GO_David_mat{"$a[1]\t$a[2]"}{$a[0]}=$a[5];
}

open(OUT,">$input_file.GO_David.mat")||die;
print OUT "Category\tTerm";

foreach my $group (@group_name){
	print OUT "\t$group";
}
print OUT "\n";

foreach my $go (sort keys %GO_David_mat){
	print OUT "$go";
	foreach my $group (@group_name){
		if(exists $GO_David_mat{$go}{$group}){
			print OUT "\t".(-log2($GO_David_mat{$go}{$group}));
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close OUT;



print "Done\n";

sub round{
	my ($s,$c)=@_;#$c, how many decimals want to keep
	my $n1=10**$c;
	#print "$c\t$n1\n";
	return $s=int(($s*$n1)+0.5)/$n1;
	
}


sub log2 {
  my $n = shift;
  return log($n)/log(2);
}




###############  End Time  ########################################################
my $Time_End="";
my $TE=time();
my $UseTime=$TE-$TS;
my $danwei="Seconds";
my $UT="";
if($UseTime>120 and $UseTime<7200){
    my $m=int($UseTime/60);
    my $s=$UseTime%60;
    $UT= "Used Time: $m minutes $s seconds";
}
elsif($UseTime<=120){
    $UT= "Used Time: $UseTime seconds";
}
elsif($UseTime>=7200){
       my $h=int($UseTime/3600);
       my $m=int(($UseTime%3600)/60);
       my $s= ($UseTime%3600)%60 ;
       $UT= "Used Time: $h hours $m minutes $s seconds";
}
 $Time_End = sub_format_datetime(localtime(time()));
 print "\nEnd Time :[$Time_End]\nTime-consuming: $UT\n";

###############################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}















