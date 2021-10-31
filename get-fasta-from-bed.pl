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

our $g_opts;
sub parse_opts{
   my $result = GetOptions(
                   "bed=s" => \$g_opts->{'bed'},#string
                   "genome=s" => \$g_opts->{'genome'},#string


                   "quiet"     => sub { $g_opts->{'verbose'} = 0 },
                   "help|?"    => \$g_opts->{'help'}
                 );
   if(!($g_opts->{'bed'}) or !($g_opts->{'genome'})){
       &pod2usage( -verbose => 1);#exit status will be 1
   }
   if($g_opts->{'help'}){
       &pod2usage( -verbose => 1);#exit status will be 1
   }
}

&parse_opts();

my $bed=$g_opts->{'bed'};
my $genome_file=$g_opts->{'genome'};

my $chr_seq=read_chrom_fasta($genome_file);

open(IN,"$bed")||die;
open(OUT,">$bed.fasta")||die;
$/="\n";
my %hash=();
my %hash2=();
my $i=0;
while(my $t=<IN>){
	chomp $t;
	$i++;
	$hash{$i}{$t}=0 if(!exists $hash2{$t});
	$hash2{$t}=0;
	#print "$t\n";
}
close IN;



foreach my $ii (sort {$a<=>$b} keys %hash){
	foreach my $t (keys %{$hash{$ii}}){
		my @a=split(/\t/,$t);
		#print "@a\n";
		print OUT ">$a[0]:$a[1]-$a[2]";
		print OUT " $a[3]" if(exists $a[3]);
		print OUT " $a[4]" if(exists $a[4]);
		print OUT "\n";
		#print $$chr_seq{$a[0]}."\n";
		#print "$a[1]\t$a[2]\n";
		my $subseq=substr($$chr_seq{$a[0]},($a[1]-1),($a[2]-$a[1])+1);
		if($a[1] eq "" or $a[2] eq "" or $a[0] eq ""){
			print "$t\n";
		}
		if(exists $a[3]){
			if($a[3] eq "-"){
				$subseq=substr($$chr_seq{$a[0]},($a[1]-1),($a[2]-$a[1])+1);
				print OUT rev($subseq)."\n";
			}
			else{
				print OUT "$subseq\n";
			}
		}
		else{
			print OUT "$subseq\n";
		}
	}
	
}

sub rev{
    my ($s)=@_;
    $s=reverse($s);

    my @a=split(//,$s);
    my $r="";
    for(my $i=0;$i<@a;$i++){

        $r.=revDNA($a[$i]);
    }
    return $r;
}



sub revDNA{
    my ($s)=@_;
    if($s ne "A" and $s ne "G" and $s ne "C" and $s ne "T" and $s ne "N" and $s ne "a" and $s ne "g" and $s ne "c" and $s ne "t" and $s ne "n" ){
        #die "Error: Wrong sequence chart $s\n";
        return "N";
        
    }
    else{
        if($s eq "A"){
            return "T";
        }
        if($s eq "T"){
            return "A";
        }
        if($s eq "C"){
            return "G";
        }
        if($s eq "G"){
            return "C";
        }
        if($s eq "N"){
            return "N";
        }
        if($s eq "a"){
            return "t";
        }
        if($s eq "t"){
            return "a";
        }
        if($s eq "c"){
            return "g";
        }
        if($s eq "g"){
            return "c";
        }
        if($s eq "n"){
            return "n";
        }
    }
}



$/=">";
sub read_chrom_fasta{
	my ($genome)=@_;
	open(IN,"$genome")||die;
	my %chr_seq=();
	$/=">";
	<IN>;
	while(my $t=<IN>){
		chomp $t;
		my @a=split(/\n/,$t);
		my $str=shift @a;
		my @aa=split(/ /,$str);
		my $subseq=join("",@a);
		#print "$subseq\n";
		$chr_seq{$aa[0]}=$subseq;
		
	}
	$/="\n";
	close IN;
	return \%chr_seq;
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




__END__

=head1 NAME

script - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

script [options] [args ...]

Options: 

   -bed           bed file 
   -genome        genome fasta file 


   -help            brief help message

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut
