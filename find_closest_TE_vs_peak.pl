use strict;
use warnings;
use File::Basename qw(basename dirname);

my $gene_file=$ARGV[0];
my $peak_file=$ARGV[1];


#open(IN,"TE_loc.txt")||die;
open(IN,$gene_file)||die;
my %TE=();
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	$TE{$a[0]}{$a[4]}="$a[1]\t$a[2]";
}

#open(IN,"h3k9me3_peak.txt")||die;
open(IN,$peak_file)||die;
my %peak=();
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	$peak{$a[0]}{$a[4]}="$a[1]\t$a[2]";
}
# run as: perl  gene_file peak_file
$peak_file=basename($peak_file);
$gene_file=~s/.txt//;
$peak_file=~s/.txt//;



open(OUT,">$gene_file.$peak_file.close_peak.txt")||die;
print OUT "chr\tstart\tend\tTE\tup_peak\tup_dist\toverlap\tdn_peak\tdn_dist\n";
foreach my $c (sort keys %TE){
	foreach my $t (sort keys %{$TE{$c}}){
		my @a_te=split(/\t/,$TE{$c}{$t});
		my $mid_te=int(($a_te[0]+$a_te[1])/2);
		my $up="-";
		my $dn="-";
		my $overlap="";
		my $up_dist=10000000000;
		my $dn_dist=10000000000;
		foreach my $p (sort keys %{$peak{$c}}){
			my @a_p=split(/\t/,$peak{$c}{$p});
			my $mid_p=int(($a_p[0]+$a_p[1])/2);
			
			if($a_p[0]<=$a_te[1] and $a_p[1]>=$a_te[0]){
				$overlap.="$p,"
			}
			elsif($a_p[1]<$a_te[0]){
				if($up_dist>$a_te[0]-$a_p[1]){
					$up_dist=$a_te[0]-$a_p[1];
					$up=$p;
				}
			}
			elsif($a_p[0]>$a_te[1]){
				if($dn_dist>$a_p[0]-$a_te[1]){
					$dn_dist=$a_p[0]-$a_te[1];
					$dn=$p;
				}
			}
			else{
				print "$p\n";
			}
		}
		$overlap=~s/,$//;
		$overlap="-" if($overlap eq "");
		$up_dist="-" if($up_dist==10000000000);
		$dn_dist="-" if($dn_dist==10000000000);
		
		print OUT "$c\t$a_te[0]\t$a_te[1]\t$t\t$up\t$up_dist\t$overlap\t$dn\t$dn_dist\n";
	}
}


# 
# open(OUT,">h3k9me3_peak.with_closest_TE.txt")||die;
# print OUT "chr\tstart\tend\tTE\tup_TE\tup_dist\tdn_TE\tdn_dist\n";
# foreach my $c (sort keys %peak){
# 	foreach my $t (sort keys %{$peak{$c}}){
# 		my @a_te=split(/\t/,$peak{$c}{$t});
# 		my $mid_te=int(($a_te[0]+$a_te[1])/2);
# 		my $up="-";
# 		my $dn="-";
# 		my $up_dist=10000000000;
# 		my $dn_dist=10000000000;
# 		foreach my $p (sort keys %{$TE{$c}}){
# 			my @a_p=split(/\t/,$TE{$c}{$p});
# 			my $mid_p=int(($a_p[0]+$a_p[1])/2);
# 			
# 			if($mid_p-$mid_te>0){
# 				if($up_dist>$mid_p-$mid_te){
# 					$up_dist=$mid_p-$mid_te;
# 					$up=$p;
# 				}
# 			}
# 			else{
# 				if($dn_dist>$mid_te-$mid_p){
# 					$dn_dist=$mid_te-$mid_p;
# 					$dn=$p;
# 				}
# 			}
# 		}
# 		$up_dist="-" if($up_dist==10000000000);
# 		$dn_dist="-" if($dn_dist==10000000000);
# 		print OUT "$c\t$a_te[0]\t$a_te[1]\t$t\t$up\t$up_dist\t$dn\t$dn_dist\n";
# 	}
# }





