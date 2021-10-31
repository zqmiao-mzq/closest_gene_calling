use strict;
use warnings;
use LWP::Simple;
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
                    "gff=s" => \$g_opts->{'gff'},#string
                    "genome=s" => \$g_opts->{'genome'},#string
                    "input=s" => \$g_opts->{'input'},#string
                    "column=i" => \$g_opts->{'column'},
                    "ORF=s" => \$g_opts->{'ORF'},#integer
                    "CDS" => \$g_opts->{'CDS'},#integer
                    "Promoter=s" => \$g_opts->{'Promoter'},#integer
                    "Tail=s" => \$g_opts->{'Tail'},#integer
                    "output=s" => \$g_opts->{'output'},
                    "fasta_format" => \$g_opts->{'fasta_format'},
                    "detail_info" => \$g_opts->{'detail_info'},

                    "quiet"     => sub { $g_opts->{'verbose'} = 0 },
                    "help|?"    => \$g_opts->{'help'}
                  );
    if(!($g_opts->{'gff'}) or !($g_opts->{'input'}) or !($g_opts->{'genome'})){
        &pod2usage( -verbose => 1);#exit status will be 1
    }
    if($g_opts->{'help'}){
        &pod2usage( -verbose => 1);#exit status will be 1
    }
}

&parse_opts();

my $gff=$g_opts->{'gff'};
my $list=$g_opts->{'input'};
my $genome_file=$g_opts->{'genome'};
my $column=0;
$column=$g_opts->{'column'}-1 if(($g_opts->{'column'}));
my $ORF="F";
my $CDS="F";
my $promoter="F";
my $tail="F";
my $out=$list;
$out= $g_opts->{'output'} if($g_opts->{'output'});

$CDS="T" if(($g_opts->{'CDS'}));

$ORF=$g_opts->{'ORF'} if(($g_opts->{'ORF'}));
$promoter=$g_opts->{'Promoter'} if(($g_opts->{'Promoter'}));
$tail=$g_opts->{'Tail'} if($g_opts->{'Tail'});



print "$gff\n$list\n$genome_file\n$promoter\n$tail\n$ORF\n$CDS\n";


my $upstream_length=$promoter;
my $downstream_length=$tail;
my %chr_name=();
my $chr_seq=read_chrom_fasta($genome_file);
my $genome=read_gff($gff,\%chr_name);

open(IN,"$list")||die;
open(PRO,">$out.CDS-$CDS.promoter-$promoter.seq")||die;
open(TAIL,">$out.CDS-$CDS.tail-$tail.seq")||die;
open(ORF,">$out.CDS-$CDS.ORF-$ORF.seq")||die;
$/="\n";
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	my $g=$a[$column];
	#print "$g\n";
	if(!exists $$genome{"gene"}{$g}){
		print "error: please check $g\n";
		next;
	}
	$g=$$genome{"gene"}{$g};
	my @info=split(/\t/,$$genome{"gene_info"}{$g});
	if($CDS eq "T"){
		$info[1]=$$genome{"gene_cds_start"}{$g} if(exists $$genome{"gene_cds_start"}{$g});
		$info[2]=$$genome{"gene_cds_end"}{$g} if(exists $$genome{"gene_cds_end"}{$g});
	}
	#print "@info\n";
	#print "$cdss\t$cdse\n";
	if($info[3] eq "+"){
		if($promoter ne "F"){
			my ($p1,$p2)=split(/,/,$promoter);
			#$p1=0-$p1;
			next if($info[1]+$p1-1<0);
			my $subseq=substr($$chr_seq{$info[0]},$info[1]+$p1-1,abs($p2-$p1));
			if($g_opts->{'detail_info'}){
				print PRO ">$g $info[0]:".($info[1]+$p1-1)."-".($info[1]+$p2-1);
				print PRO "|info=$info[0]:$info[1]-$info[2]|$info[3]";
				print PRO "|desc=".$$genome{"gene_desc"}{$g} if(exists $$genome{"gene_desc"}{$g});
				print PRO "|common_name=".$$genome{"gene2common"}{$g} if(exists $$genome{"gene2common"}{$g});
				print PRO "\n";
			}
			else{
				print PRO ">$g\n";
			}
			if($g_opts->{'fasta_format'}){
				$subseq=fasta_format($subseq);
			}
			$subseq=~s/\n$//;
			print PRO "$subseq\n";
		}
		if($tail ne "F"){
			my ($t1,$t2)=split(/,/,$tail);
			#$t1=0-$t1;
			next if($info[2]+1+$t1<0);
			my $subseq=substr($$chr_seq{$info[0]},$info[2]+1+$t1,abs($t2-$t1));
			if($g_opts->{'detail_info'}){
				print TAIL ">$g $info[0]:".($info[2]+1+$t1)."-".($info[2]+1+$t2);
				print TAIL "|$info[0]:$info[1]-$info[2]|$info[3]";
				print TAIL "|".$$genome{"gene_desc"}{$g} if(exists $$genome{"gene_desc"}{$g});
				print TAIL "\n";
			}
			else{
				print TAIL ">$g\n";
			}
			if($g_opts->{'fasta_format'}){
				$subseq=fasta_format($subseq);
			}
			$subseq=~s/\n$//;
			print TAIL "$subseq\n";
		}
		if($ORF ne "F"){
			my $subseq="";
			print "$ORF------\n";
			if($ORF eq "All"){
				my $ORF_l=$info[2]-$info[1]+1;
				next if($info[1]-1<0);
				$subseq=substr($$chr_seq{$info[0]},$info[1]-1,$ORF_l); 
				print "$ORF_l=$info[2]-$info[1]\n$subseq\n\n";
			}
			elsif(!($ORF=~/,/)){
				next if($info[1]-1<0);
				$subseq=substr($$chr_seq{$info[0]},$info[1]-1,$ORF); 
			}
			else{
				my ($o1,$o2)=split(/,/,$ORF);
				next if($info[1]-1+$o1<0);
				$subseq=substr($$chr_seq{$info[0]},$info[1]-1+$o1,abs($o2+$info[2]-($o1+$info[1])));
			}
			
			if($g_opts->{'detail_info'}){
				print ORF ">$g";
				print ORF "|$info[0]:$info[1]-$info[2]|$info[3]";
				print ORF "|".$$genome{"gene_desc"}{$g} if(exists $$genome{"gene_desc"}{$g});
				print ORF "\n";
			}
			else{
				print ORF ">$g\n";
			}
			if($g_opts->{'fasta_format'}){
				$subseq=fasta_format($subseq);
			}
			$subseq=~s/\n$//;
			print ORF "$subseq\n";
		}
		
	}
	else{
		if($promoter ne "F"){
			my ($p1,$p2)=split(/,/,$promoter);
			#$p1=0-$p1;
			next if($info[2]-$p2<0);
			my $subseq=substr($$chr_seq{$info[0]},$info[2]-$p2,abs($p2-$p1));
			$subseq=rev($subseq);
			if($g_opts->{'detail_info'}){
				print PRO ">$g $info[0]:".($info[1]+$p1-1)."-".($info[1]+$p2-1);
				print PRO "|$info[0]:$info[1]-$info[2]|$info[3]";
				print PRO "|".$$genome{"gene_desc"}{$g} if(exists $$genome{"gene_desc"}{$g});
				print PRO "\n";
			}
			else{
				print PRO ">$g\n";
			}
			if($g_opts->{'fasta_format'}){
				$subseq=fasta_format($subseq);
			}
			$subseq=~s/\n$//;
			print PRO "$subseq\n";
		}
		if($tail ne "F"){
			my ($t1,$t2)=split(/,/,$tail);
			#$t1=0-$t1;
			next if($info[1]-$t2-1<0);
			my $subseq=substr($$chr_seq{$info[0]},$info[1]-$t2-1,abs($t2-$t1));
			$subseq=rev($subseq);
			if($g_opts->{'detail_info'}){
				print TAIL ">$g $info[0]:".($info[1]-$t2-1)."-".($info[1]-$t1-1);
				print TAIL "|$info[0]:$info[1]-$info[2]|$info[3]";
				print TAIL "|".$$genome{"gene_desc"}{$g} if(exists $$genome{"gene_desc"}{$g});
				print TAIL "\n";
			}
			else{
				print TAIL ">$g\n";
			}
			if($g_opts->{'fasta_format'}){
				$subseq=fasta_format($subseq);
			}
			$subseq=~s/\n$//;
			print TAIL "$subseq\n";
		}
		if($ORF ne "F"){
			my $subseq="";
			if($ORF eq "All"){
				my $ORF_l=$info[2]-$info[1]+1;
				next if($info[1]-1<0);
				$subseq=substr($$chr_seq{$info[0]},$info[1]-1,$ORF_l); 
				print "$ORF_l=$info[2]-$info[1]\n$subseq\n\n";
			}
			elsif(!($ORF=~/,/)){
				next if($info[2]-$ORF<0);
				$subseq=substr($$chr_seq{$info[0]},$info[2]-$ORF,$ORF); 
			}
			else{
				my ($o1,$o2)=split(/,/,$ORF);
				next if($info[1]-$o2<0);
				$subseq=substr($$chr_seq{$info[0]},$info[1]-$o2,abs($info[2]-$o1-$info[1]+$o2));
			}
			
			$subseq=rev($subseq);
			if($g_opts->{'detail_info'}){
				print ORF ">$g";
				print ORF "|$info[0]:$info[1]-$info[2]|$info[3]";
				print ORF "|".$$genome{"gene_desc"}{$g} if(exists $$genome{"gene_desc"}{$g});
				print ORF "\n";
			}
			else{
				print ORF ">$g\n";
			}
			if($g_opts->{'fasta_format'}){
				$subseq=fasta_format($subseq);
			}
			$subseq=~s/\n$//;
			print ORF "$subseq\n";
		}
	
	
	
	}
	
}

`rm $out.CDS-$CDS.promoter-$promoter.seq` if($promoter eq "F");
`rm $out.CDS-$CDS.ORF-$ORF.seq` if($ORF eq "F");
`rm $out.CDS-$CDS.tail-$tail.seq` if($tail eq "F");





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

sub fasta_format{
    my ($s)=@_;
	my $n=int((length($s))/60);
	my $res="";
	for(my $i=0;$i<=$n;$i++){
		my $tmp=substr($s,$i*60,60);
		$res.=$tmp."\n" if($tmp ne "");
	}
	return $res;
	
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


sub read_gff{
    my ($gff,$chr_name)=@_;
    open(IN,"$gff")||die "can't open $gff\n";
    my %chr_length=();
	$/="\n";
    #my %chr_name=();
    my %gene_cds_start=();
    my %gene_cds_end=();
    my %gene_info=();
    my %gene=();
    my %gene_name2gene_id=();
    my %mRNA=();
    my %gene2ID=();
    my %ID2gene=();
    my %gene2common=();
	my %exon=();
	my %CDS=();
	my %promoter_region=();
	my %tail_region=();
    my %all_gene_name=();
    my %start2gene=();
    my %end2gene=();
    my %cds_start2gene=();
    my %cds_end2gene=();
    my %gene_cds=();
    my %gene_desc=();
    my %gene_hash=();
    my %exclude=();
    while(my $t=<IN>){
            chomp $t;
            next if($t=~/^\#/);
            my @a=split(/\t/,$t);
            if($a[2] eq "region" or $a[2] eq "chromosome"){
				if(exists $$chr_name{$a[0]}){
					$a[0]=$$chr_name{$a[0]};
				}
                if(exists $chr_length{$a[0]}){
                    $chr_length{$a[0]}=$a[4] if($a[4]>$chr_length{$a[0]});
                }
                else{
                     $chr_length{$a[0]}=$a[4];
                }
                #print "$a[0]\t$a[4]\n";
            }
        	if($a[2] eq "ncRNA" or $a[2] eq "tRNA" or $a[2] eq "rRNA"){
        		$a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }    
                $exclude{$tmphash{"ID"}}=0 if(exists $tmphash{"ID"});
                $exclude{$tmphash{"Parent"}}=0 if(exists $tmphash{"Parent"});
        	}
        	
        	
        
            if($a[2] eq "CDS"){
                
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }      

                
                my $mrna_ID=$tmphash{"Parent"};

            
                if(exists $mRNA{$mrna_ID}){
                    next if(!exists $ID2gene{$mRNA{$mrna_ID}});
                    $CDS{$a[0]}{$ID2gene{$mRNA{$mrna_ID}}}{$a[3]."\t".$a[4]}=0;
                    if(exists $gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}}){ 
                        $gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}}=$a[3] if($a[3]<$gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}});
                        $gene_cds_end{$ID2gene{$mRNA{$mrna_ID}}}=$a[4] if($a[4]>$gene_cds_end{$ID2gene{$mRNA{$mrna_ID}}});
                        #print "$mRNA{$mrna_ID}\t$a[3]\n";
                        
                    }
                    else{
                    	
                        $gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}}=$a[3];
                        $gene_cds_end{$ID2gene{$mRNA{$mrna_ID}}}=$a[4];
                        #print "$mrna_ID\t$a[3]\n";
                    }
                }
                else{
                	next if(!exists $ID2gene{$mrna_ID});
                	$CDS{$a[0]}{$ID2gene{$mrna_ID}}{$a[3]."\t".$a[4]}=0;
                    if(exists $gene_cds_start{$ID2gene{$mrna_ID}}){
                        
                        $gene_cds_start{$ID2gene{$mrna_ID}}=$a[3] if($a[3]<$gene_cds_start{$ID2gene{$mrna_ID}});
                        $gene_cds_end{$ID2gene{$mrna_ID}}=$a[4] if($a[4]>$gene_cds_end{$ID2gene{$mrna_ID}});
                        #print "$mRNA{$mrna_ID}\t$a[3]\n";
                    }
                    else{
                        
                        $gene_cds_start{$ID2gene{$mrna_ID}}=$a[3];
                        $gene_cds_end{$ID2gene{$mrna_ID}}=$a[4];
                        #print "$mrna_ID\t$a[3]\n";
                    }
                }
            
            
            
        
            }
            if($a[2] eq "mRNA"){
                
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }      
                
                my $ID=$tmphash{"Parent"};
                my $mrna_ID=$tmphash{"ID"};

                $mRNA{$mrna_ID}=$ID;

        
            }
            
            if($a[2] eq "exon"){
                
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }
                
                
                my $mrna_ID=$tmphash{"Parent"};

                
				if(exists $$chr_name{$a[0]}){
					$a[0]=$$chr_name{$a[0]};
				}

                
                if(exists $mRNA{$mrna_ID}){
                	next if(!exists $ID2gene{$mRNA{$mrna_ID}});
                    $exon{$a[0]}{$ID2gene{$mRNA{$mrna_ID}}}{$a[3]."\t".$a[4]}=0;
                }
                else{
                	next if(!exists $ID2gene{$mrna_ID});
                    $exon{$a[0]}{$ID2gene{$mrna_ID}}{$a[3]."\t".$a[4]}=0;
                }               
        
            }

            if($a[2] eq "gene" or $a[2] eq "pseudogene"){
                if($chr_name){
                    if(exists $$chr_name{$a[0]}){
                        $a[0]=$$chr_name{$a[0]};
                    }
                }
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }
                if($a[8]=~/GeneID:(\d+)/){
                	$tmphash{"entrz_id"}=$1; 
                }
                if($a[8]=~/HGNC:(\d+)/){
                	$tmphash{"hgnc_id"}=$1; 
                }
                
                if($a[0] eq "ChrM" or $a[0] eq "NC_012920.1" or $a[0] eq "mito_A_nidulans_FGSC_A4" or $a[0] eq "Ca19-mtDNA" or $a[0] eq "NC_010339.1" or $a[0] eq "NC_005089.1"){
                	$exclude{$tmphash{"ID"}}=0 if(exists $tmphash{"ID"});
                	$exclude{$tmphash{"Parent"}}=0 if(exists $tmphash{"Parent"});
                	
                }
                if(exists $tmphash{"pseudo"} and $tmphash{"pseudo"} eq "true"){
                	$exclude{$tmphash{"ID"}}=0 if(exists $tmphash{"ID"});
                	$exclude{$tmphash{"Parent"}}=0 if(exists $tmphash{"Parent"});
                }
              	
              	
            
                my $ID=$tmphash{"ID"};
                my $gene_id=$tmphash{"entrz_id"} if(exists $tmphash{"entrz_id"});
            
                my $gene_name=$ID;
                $gene_name=$tmphash{"Name"} if(exists $tmphash{"Name"});
                $gene_name=$tmphash{"gene_name"} if(exists $tmphash{"gene_name"});


                my $common="";
                $common=$tmphash{"gene"} if(exists $tmphash{"gene"});
                $common=$tmphash{"Gene"} if(exists $tmphash{"Gene"});

                $gene{$gene_name}=$gene_name;
                

                $gene_name2gene_id{$gene_name}=$gene_id if(defined $gene_id);
                
                $all_gene_name{$gene_name}=0;
                
                $gene2ID{$gene_id}=$ID if(defined $gene_id);
                
                $gene2ID{$gene_name}=$ID;
                
                $ID2gene{$ID}=$gene_name;
                
                $gene2common{$gene_id}=$common if($common ne "" and defined $gene_id);
                
                $gene{$gene_id}=$gene_name if(defined $gene_id);
                
                $gene2common{$gene_name}=$common if($common ne "");
                $gene{$common}=$gene_name if($common ne "");
                $gene{uc($common)}=$gene_name if($common ne "");

                %{$gene_hash{$gene_name}}=%tmphash;
                
                $gene_info{$gene_name}="$a[0]\t$a[3]\t$a[4]\t$a[6]";
                $gene_info{$gene_id}="$a[0]\t$a[3]\t$a[4]\t$a[6]" if(defined $gene_id);
                if($a[6] eq "+"){
                	$start2gene{$a[0]}{$a[6]}{$a[3]}{$gene_name}=0;
                	$end2gene{$a[0]}{$a[6]}{$a[4]}{$gene_name}=0;
                }
                else{
                	$start2gene{$a[0]}{$a[6]}{$a[4]}{$gene_name}=0;
                	$end2gene{$a[0]}{$a[6]}{$a[3]}{$gene_name}=0;
                }
                my $alias=$tmphash{"Alias"} if(exists $tmphash{"Alias"});
                my $note=$tmphash{"Note"} if(exists $tmphash{"Note"});
                my $desc=$tmphash{"description"} if(exists $tmphash{"description"});

                
                if(defined $alias and $alias ne ""){
                    my @aa=split(/,/,$alias);
                    foreach my $tt (@aa){
                            if(!exists $gene{$tt}){
                                    $gene{$tt}=$gene_name;
                            }
                    }
                }
                
                
                if(defined $note){
                    $gene_desc{$gene_name}="$note";
                    $gene_desc{$gene_id}="$note" if(defined $gene_id);
                }
                
                if(defined $desc){
                    $gene_desc{$gene_name}="$desc";
                    $gene_desc{$gene_id}="$desc" if(defined $gene_id);
                }
            }
        
    }
    
    
    foreach my $g (keys %all_gene_name){
    	
    	my @info=split(/\t/,$gene_info{$g});
    	if(exists $gene_cds_start{$g}){
    		$gene_cds{$info[0]}{$g}=$gene_cds_start{$g}."\t".$gene_cds_end{$g};
    		if($info[3] eq "+"){
    			$cds_start2gene{$info[0]}{$info[3]}{$gene_cds_start{$g}}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$gene_cds_end{$g}}{$g}=0;
    		}
    		else{
    			$cds_start2gene{$info[0]}{$info[3]}{$gene_cds_end{$g}}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$gene_cds_start{$g}}{$g}=0;
    		}
    		
    	}
    	else{
    		$gene_cds{$info[0]}{$g}=$info[1]."\t".$info[2];
			if($info[3] eq "+"){
    			$cds_start2gene{$info[0]}{$info[3]}{$info[1]}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$info[2]}{$g}=0;	
    		}
    		else{
    			$cds_start2gene{$info[0]}{$info[3]}{$info[2]}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$info[1]}{$g}=0;
    		}
    	}
    }
    
    

    my %genome=();
    %{$genome{"gene"}}=%gene;
    %{$genome{"exon"}}=%exon;
    %{$genome{"gene_cds"}}=%gene_cds;
    
    %{$genome{"promoter_region"}}=%promoter_region;
    %{$genome{"tail_region"}}=%tail_region;
    %{$genome{"exon"}}=%exon;
    %{$genome{"CDS"}}=%CDS;
    %{$genome{"ID2gene"}}=%ID2gene;
    %{$genome{"all_gene_name"}}=%all_gene_name;
    
    %{$genome{"chr_length"}}=%chr_length;
    %{$genome{"gene_info"}}=%gene_info;
    %{$genome{"gene_cds_start"}}=%gene_cds_start;
    %{$genome{"gene_cds_end"}}=%gene_cds_end;
    
    %{$genome{"cds_start2gene"}}=%cds_start2gene;
    %{$genome{"cds_end2gene"}}=%cds_end2gene;
    %{$genome{"start2gene"}}=%start2gene;
    %{$genome{"end2gene"}}=%end2gene;
    
    %{$genome{"gene2common"}}=%gene2common;
    %{$genome{"gene2ID"}}=%gene2ID;
    %{$genome{"gene_name2gene_id"}}=%gene_name2gene_id;
    %{$genome{"gene_hash"}}=%gene_hash;
    %{$genome{"gene_desc"}}=%gene_desc;
    %{$genome{"exclude"}}=%exclude;

    
    return \%genome;


}



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
		$chr_seq{$aa[0]}=$subseq;
	}
	$/="\n";
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






# print("\n$optionX\n");
# print($optionY); #if not defined in command line, it will be undefined
# print($g_opts->{"optionB"});
# 
# foreach my $key (keys %{$g_opts}){
#   if(!$g_opts->{$key}) {next;} 
#   print($key . "=" . $g_opts->{$key} . "\n");
#   }



__END__

=head1 NAME

script - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

script [options] [args ...]

Options: 

   -gff           gff file 
   -genome        genome fasta file 
   
   -input         ur input list
   -column        int, which column u want to handle, count from 1, defined is 1
   
   -ORF           T or F, T means get the ORF sequence, F means not, defined is F, All means the whole gene sequence,like -1000,1000 or 400
   -Promoter      the length of the region before start,like -2000,0, F means do not get the sequence of promoter, defined is F
   -Tail          the length of the region after end, like 0,1000, F means do not get the sequence of tail, defined is F
   
   -CDS           if you want to use CDS position
   -fasta_format  if you want to output fasta_format
   -detail_info   if you want detail gene information
   
   -output        output name

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
