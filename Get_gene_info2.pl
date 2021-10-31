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



my %gff=();

$gff{"AN"}="/Users/zqmiao/Documents/data/genome/A_nidulans/A_nidulans_FGSC_A4_version_s10-m04-r16_features.gff";
$gff{"AF"}="/Users/zqmiao/Documents/data/genome/A_fumigatus/A_fumigatus_Af293_current_features.gff";
$gff{"CA"}="/Users/zqmiao/Documents/data/genome/c_albicans/C_albicans_SC5314_A21_current_features.modified.gff";
$gff{"CG"}="/Users/zqmiao/Documents/data/genome/C_glabrata/C_glabrata_CBS138_current_features.gff";
$gff{"Mouse"}="";
$gff{"Human"}="/Users/zqmiao/Documents/data/genome/human/GCF_000001405.28_GRCh38.p2_genomic.gff";
$gff{"ZebraFish"}="";





my %NCBI_ID=();
my %NCBI_Synonyms_ID=();
my %Info=();
my %NCBI_ID2KEGG_ID=();
open(IN,"$gene_list")||die "can't read $gene_list\n";
print "$gene_list\n";
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

my @input_mat=();
while(my $t=<IN>){
	chomp $t;
	my @a=split(/\t/,$t);
	#push(@{$input_mat{$a[1]}},$a[0]);

	for(my $i=0; $i<@a;$i++){
		push(@input_mat,$a[$i]);
	}
}
close IN;
#close OUT;



my %chr_name=();
my $genome=read_gff($gff{$specis},\%chr_name);




my $path=dirname($input_file);

open(OUT,">$input_file.gene_info")||die "$input_file.gene_info";
print OUT "Input1\tTax_ID\tNCBI_ID\tSymbol\tLocusTag\tSynonyms\tKEGG_ID\tLab_ID\tCommonName\tdescription\n";


foreach my $g (@input_mat){
	if(exists $NCBI_ID{$g}){
		print OUT "$g\t$Info{$NCBI_ID{$g}}";
	}
	else{
		print OUT "$g\t-\t-\t-\t-\t-\t-\t-\t-\t-";
	}
	
	if(exists $$genome{"gene_info"}{$g}){
		print OUT "\t".$$genome{"gene_info"}{$g};
	}
	else{
		print OUT "\t-\t-\t-\t-";
	}
	if(exists $$genome{"gene_desc"}{$g}){
		print OUT "\t".$$genome{"gene_desc"}{$g}."\n";
	}
	else{
		print OUT "\t-\n";
	}
}

close OUT;


sub read_gff{
    my ($gff,$chr_name)=@_;
    open(IN,"$gff")||die "can't open $gff\n";
    my %chr_length=();

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
    my %gene_start=();
    my %gene_end=();
    my %cds_start2gene=();
    my %cds_end2gene=();
    my %gene_cds=();
    my %gene_desc=();
    my %gene_hash=();
    my %exclude=();
    my %genebody=();
    $/="\n";
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
        	if($a[2] eq "ncRNA" or $a[2] eq "tRNA" or $a[2] eq "rRNA" or $a[2] eq "snRNA" or $a[2] eq "sRNA" or $a[2] eq "snoRNA" or $a[2] eq "siRNA" or $a[2] eq "scaRNA" or $a[2] eq "piRNA" or $a[2] eq "centromere" or $a[2] eq "repeat_region" or $a[2] eq "retrotransposon" or $a[2] eq "long_terminal_repeat" or $a[2] eq "blocked_reading_frame"){
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
                $exclude{$tmphash{"Name"}}=0 if(exists $tmphash{"Name"});
                $exclude{$tmphash{"name"}}=0 if(exists $tmphash{"name"});
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
                if($a[6] eq "+"){
                	$gene_start{$gene_name}=$a[3];
                	$gene_end{$gene_name}=$a[4];
                }
                else{
                	$gene_start{$gene_name}=$a[4];
                	$gene_end{$gene_name}=$a[3];
            	}
                $gene_info{$gene_name}="$a[0]\t$a[3]\t$a[4]\t$a[6]";
                
                # for(my $gb=$a[3];$gb<=$a[4];$gb++){
#                 	$genebody{$a[0]}{$gb}{$gene_name}=0;
#                 }
                
                $genebody{$a[0]}{$gene_name}=0;
                
                if(defined $gene_id){
                	$gene_info{$gene_id}="$a[0]\t$a[3]\t$a[4]\t$a[6]";
                	
                }
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
    
	my %gene_info_CDS=();
    foreach my $g (keys %all_gene_name){
    	
    	my @info=split(/\t/,$gene_info{$g});
		
    	if(exists $gene_cds_start{$g}){
    		$gene_cds{$info[0]}{$g}=$gene_cds_start{$g}."\t".$gene_cds_end{$g};
    		$gene_info_CDS{$g}="$info[0]\t$gene_cds_start{$g}\t$gene_cds_end{$g}\t$info[3]";
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
    		$gene_info_CDS{$g}=$gene_info{$g};
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
    %{$genome{"gene_info_CDS"}}=%gene_info_CDS;
    %{$genome{"gene_cds_start"}}=%gene_cds_start;
    %{$genome{"gene_cds_end"}}=%gene_cds_end;
    %{$genome{"gene_start"}}=%gene_start;
    %{$genome{"gene_end"}}=%gene_end;
    
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
    %{$genome{"genebody"}}=%genebody;

    
    return \%genome;


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















