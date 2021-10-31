




sub split_arr{
	my ($a,$threads)=@_;
	my $n=@$a;
	print "$n\n";
	my $sub_n=int($n/$threads)+1;
	my @res=();
	my $j=0;
	for(my $i=0;$i<$n;$i++){
		push(@{$res[$j]},$$a[$i]);
		$j++ if(($i+1)%$sub_n==0);
	}
	return \@res;
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
    if($s ne "A" and $s ne "G" and $s ne "C" and $s ne "T" and $s ne "N" ){
        die "Error: Wrong sequence chart $s\n";
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

    }

}




sub score_pwm{
	my ($seq,$pwm)=@_;
 	my @a=split(//,$seq);
 	my %res=();

	for(my $i=0;$i<@a;$i++){
		#print "$i\n" if($i%10000==0);
		foreach my $l (keys %{$$pwm{"length"}}){
			next if($i+$l>@a);
			foreach my $m (keys %{$$pwm{"length"}{$l}}){
				my $score1=0;
				my $score2=0;
				my $subseq1="";
				my $subseq2="";
				for(my $j=0;$j<$l;$j++){
					my $c=$a[$i+$j];	
					if(!defined $$pwm{"matrix"}{$m}{$c}[$j]){
						$score1="";
						last;
					}
 					print "i-$i\tl-$l\tm-$m\tc-$c\tj-$j\t".$$pwm{"matrix"}{$m}{$c}[$j]."\n" if(!defined $$pwm{"matrix"}{$m}{$c}[$j]);

					$score1+=$$pwm{"matrix"}{$m}{$c}[$j];
					$subseq1.=$c;
					
					my $c2=revDNA($a[$i+$l-1-$j]);
					if(!defined $$pwm{"matrix"}{$m}{$c2}[$j]){
						$score2="";
						last;
					}
					$score2+=$$pwm{"matrix"}{$m}{$c2}[$j];
					$subseq2.=$c2;
				}
				#my $str=$i."\t".($i+$l-1)."\t".$subseq1."\t+\t".round(($score1-$$pwm{"min"}{$m})/($$pwm{"max"}{$m}-$$pwm{"min"}{$m}),4)."\t$m";
				#print $str."\n";
				next if($score1 eq "" or $score2 eq "");
				if(round(($score1-$$pwm{"min"}{$m})/($$pwm{"max"}{$m}-$$pwm{"min"}{$m}),4)>=0.9){
					my $str=$i."\t".($i+$l-1)."\t".$subseq1."\t+\t".round(($score1-$$pwm{"min"}{$m})/($$pwm{"max"}{$m}-$$pwm{"min"}{$m}),4)."\t$m";
					$res{$str}=0;
					#print "$str\n";
				}

				if(round(($score2-$$pwm{"min"}{$m})/($$pwm{"max"}{$m}-$$pwm{"min"}{$m}),4)>=0.9){
					my $str=$i."\t".($i+$l-1)."\t".$subseq2."\t-\t".round(($score2-$$pwm{"min"}{$m})/($$pwm{"max"}{$m}-$$pwm{"min"}{$m}),4)."\t$m";
					$res{$str}=0;
					#print "$str\n";
				}

			}
		}
	}
	return \%res;
}



sub min{
	my ($a)=@_;
	my $min=$$a[0];
	foreach my $s (@$a){
		if($min>$s){
			$min=$s;
		}
	}
	return $min;
}

sub max{
	my ($a)=@_;
	my $max=$$a[0];
	#print "$max=1\n";
	foreach my $s (@$a){
		if($max<$s){
		    #print "max:\t$max\t$s\n";
			$max=$s;
		}
	}
	return $max;
}

sub round{
	my ($s1,$c)=@_;#$c, how many decimals want to keep
	my $n1=10**$c;
	my $s2=abs($s1);
	#print "\n$s1\t$s2\t$c\t$n1\n";
	return int($s1*$n1+0.5)/$n1;
	
}


sub read_pwm{
	my ($pwm_file)=@_;
	open(IN,"$pwm_file")||die;
	my %pwm=();
	my %bin_hash=();
	my %pwm_min=();
	my %pwm_max=();
	$/=">";
	<IN>;
	while(my $t=<IN>){
		chomp $t;
		my @a=split(/\n/,$t);
		my $na=shift @a;
		
		
		my %chart=();
		@{$chart{"A"}}=split(/\t/,$a[0]);
		@{$chart{"C"}}=split(/\t/,$a[1]);
		@{$chart{"G"}}=split(/\t/,$a[2]);
		@{$chart{"T"}}=split(/\t/,$a[3]);
		%{$pwm{"matrix"}{$na}}=%chart;
		my $length=@{$chart{"A"}};
		$pwm{"length"}{$length}{$na}=0;
		#print "@{$chart{\"A\"}}\n";
		my $s_min=0;
		my $s_max=0;
		for(my $i=0;$i<$length;$i++){
			my @tmp=();
			foreach my $c (keys %chart){
				#print "$c\t$chart{$c}[$i]\n";
				push(@tmp,$chart{$c}[$i]);

			}
			$s_min+=min(\@tmp);
			$s_max+=max(\@tmp);
		}
		$pwm{"min"}{$na}=$s_min;
		$pwm{"max"}{$na}=$s_max;
		
	}
	return \%pwm;	
	
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
		$chr_seq{$aa[0]}=uc($subseq);
	}
	return \%chr_seq;
}



sub read_gff{

#     %{$genome{"gene"}}=%gene; #$$genome{"gene"}{all_kind_gene_ID/name}=genename
#     %{$genome{"gene_cds"}}=%gene_cds; #$$genome{"gene_cds"}{gene_ID/genename}=chr\tcds_start\tcds_end
#     
#     %{$genome{"exon"}}=%exon; #$$genome{"exon"}{chr}{gene}{start-end}=0
#     %{$genome{"CDS"}}=%CDS; #$$genome{"CDS"}{chr}{gene}{start\tend}=0
#     %{$genome{"ID2gene"}}=%ID2gene; #$$genome{"ID2gene"}{gene_id}=genename
#     %{$genome{"all_gene_name"}}=%all_gene_name; #$$genome{"all_gene_name"}{genename}=0
#     
#     %{$genome{"chr_length"}}=%chr_length;# $$genome{"chr_length"}{chr}=length
#     %{$genome{"gene_info"}}=%gene_info; #$$genome{"gene_info"}{gene_ID/genename}=chr\tstart\tend\tstrand
#     %{$genome{"gene_info2"}}=%gene_info2; #$$genome{"gene_info2"}{chr}{gene_ID/genename}=chr\tstart\tend\tstrand
#     %{$genome{"gene_cds_start"}}=%gene_cds_start;# $$genome{"gene_cds_start"}{gene}=start
#     %{$genome{"gene_cds_end"}}=%gene_cds_end;# $$genome{"gene_cds_start"}{gene}=end
#     
#     %{$genome{"cds_start2gene"}}=%cds_start2gene; # $$genome{"cds_start2gene"}{chr}{strand}{start}{genename}=0 #direction considered
#     %{$genome{"cds_end2gene"}}=%cds_end2gene; # $$genome{"cds_start2gene"}{chr}{strand}{end}{genename}=0 #direction considered
#     %{$genome{"start2gene"}}=%start2gene; # $$genome{"start2gene"}{chr}{strand}{end}{genename}=0 #direction considered
#     %{$genome{"end2gene"}}=%end2gene; # $$genome{"end2gene"}{chr}{strand}{end}{genename}=0 #direction considered
#     
#     %{$genome{"gene2common"}}=%gene2common; # $$genome{"gene2common"}{gene_ID/name}=common_name 
#     %{$genome{"gene2ID"}}=%gene2ID; # $$genome{"gene2ID"}{gene_ID/name}=gene_ID 
#     %{$genome{"gene_name2gene_id"}}=%gene_name2gene_id; # $$genome{"gene_name2gene_id"}{genename}=gene_ID
#     %{$genome{"gene_hash"}}=%gene_hash;  # $$genome{"gene_hash"}{genename}{all_gene_feature}=feature
#     %{$genome{"gene_desc"}}=%gene_desc;  # $$genome{"gene_desc"}{gene_ID/name}=desc/note
#     %{$genome{"exclude"}}=%exclude;  # $$genome{"exclude"}{gene_ID/name}=0


    my ($gff,$chr_name)=@_;
    open(IN,"$gff")||die "can't open $gff\n";
    my %chr_length=();
	use URI::Escape;
    #my %chr_name=();
    my %gene_cds_start=();
    my %gene_cds_end=();
    my %gene_info=();
    my %gene_info2=();
    my %gene=();
    my %gene_name2gene_id=();
    my %mRNA=();
    my %gene2ID=();
    my %ID2gene=();
    my %gene2common=();
	my %exon=();
	my %CDS=();
    my %all_gene_name=();
    my %start2gene=();
    my %end2gene=();
    my %cds_start2gene=();
    my %cds_end2gene=();
    my %gene_cds=();
    my %gene_desc=();
    my %gene_hash=();
    my %exclude=();
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
        	if($a[2] eq "ncRNA" or $a[2] eq "tRNA" or $a[2] eq "rRNA" or $a[2] eq "snRNA" or $a[2] eq "sRNA" or $a[2] eq "snoRNA" or $a[2] eq "siRNA" or $a[2] eq "scaRNA" or $a[2] eq "piRNA"){
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
              	
              	
                my $gene_name="";
                $gene_name=$tmphash{"Name"} if(exists $tmphash{"Name"});
                $gene_name=$tmphash{"gene_name"} if(exists $tmphash{"gene_name"});

                my $ID=$tmphash{"ID"};
                my $gene_id=$tmphash{"entrz_id"} if(exists $tmphash{"entrz_id"});
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
                $gene_info2{$a[0]}{$gene_name}="$a[0]\t$a[3]\t$a[4]\t$a[6]";
                $gene_info2{$a[0]}{$gene_id}="$a[0]\t$a[3]\t$a[4]\t$a[6]" if(defined $gene_id);
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
    
    close IN;
    
    foreach my $g (keys %all_gene_name){
    	
    	my @info=split(/\t/,$gene_info{$g});
    	if(exists $gene_cds_start{$g}){
    		$gene_cds{$g}=$info[0]."\t".$gene_cds_start{$g}."\t".$gene_cds_end{$g}."\t".$info[3];
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
    		$gene_cds{$g}=$info[0]."\t".$info[1]."\t".$info[2]."\t".$info[3];
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
    %{$genome{"gene"}}=%gene; #$$genome{"gene"}{all_kind_gene_ID/name}=genename
    %{$genome{"gene_cds"}}=%gene_cds; #$$genome{"gene_cds"}{gene_ID/genename}=chr\tcds_start\tcds_end
    
    %{$genome{"exon"}}=%exon; #$$genome{"exon"}{chr}{gene}{start-end}=0
    %{$genome{"CDS"}}=%CDS; #$$genome{"CDS"}{chr}{gene}{start\tend}=0
    %{$genome{"ID2gene"}}=%ID2gene; #$$genome{"ID2gene"}{gene_id}=genename
    %{$genome{"all_gene_name"}}=%all_gene_name; #$$genome{"all_gene_name"}{genename}=0
    
    %{$genome{"chr_length"}}=%chr_length;# $$genome{"chr_length"}{chr}=length
    %{$genome{"gene_info"}}=%gene_info; #$$genome{"gene_info"}{gene_ID/genename}=chr\tstart\tend\tstrand
    %{$genome{"gene_info2"}}=%gene_info2; #$$genome{"gene_info2"}{chr}{gene_ID/genename}=chr\tstart\tend\tstrand
    %{$genome{"gene_cds_start"}}=%gene_cds_start;# $$genome{"gene_cds_start"}{gene}=start
    %{$genome{"gene_cds_end"}}=%gene_cds_end;# $$genome{"gene_cds_start"}{gene}=end
    
    %{$genome{"cds_start2gene"}}=%cds_start2gene; # $$genome{"cds_start2gene"}{chr}{strand}{start}{genename}=0 #direction considered
    %{$genome{"cds_end2gene"}}=%cds_end2gene; # $$genome{"cds_start2gene"}{chr}{strand}{end}{genename}=0 #direction considered
    %{$genome{"start2gene"}}=%start2gene; # $$genome{"start2gene"}{chr}{strand}{end}{genename}=0 #direction considered
    %{$genome{"end2gene"}}=%end2gene; # $$genome{"end2gene"}{chr}{strand}{end}{genename}=0 #direction considered
    
    %{$genome{"gene2common"}}=%gene2common; # $$genome{"gene2common"}{gene_ID/name}=common_name 
    %{$genome{"gene2ID"}}=%gene2ID; # $$genome{"gene2ID"}{gene_ID/name}=gene_ID 
    %{$genome{"gene_name2gene_id"}}=%gene_name2gene_id; # $$genome{"gene_name2gene_id"}{genename}=gene_ID
    %{$genome{"gene_hash"}}=%gene_hash;  # $$genome{"gene_hash"}{genename}{all_gene_feature}=feature
    %{$genome{"gene_desc"}}=%gene_desc;  # $$genome{"gene_desc"}{gene_ID/name}=desc/note
    %{$genome{"exclude"}}=%exclude;  # $$genome{"exclude"}{gene_ID/name}=0

    
    return \%genome;


}


sub avg{
    my ($a)=@_;
    my $n=@$a;
    my $sum=0;
    foreach my $t (@$a){
        $sum+=$t;
    }
    return $sum/$n;

}


sub round{
	my ($s,$c)=@_;#$c, how many decimals want to keep
	my $n1=10**$c;
	#print "$c\t$n1\n";
	return $s=int(($s*$n1)+0.5)/$n1;
}


sub pfm2pwm{
    my ($h,$pseudocount,$acgt)=@_;
    die if($h eq "" or $pseudocount eq "" or $acgt eq "");
    
    my $sum=0;

    foreach my $tt (sort {$a<=>$b} keys %{$h}){
        $sum+=$$h{$tt}[0];
    }
    #print "$sum\n";
    my @background=split(/,/,$acgt);

    foreach my $tt (sort {$a<=>$b} keys %{$h}){
        for(my $j=0;$j<@{$$h{$tt}};$j++){
        	
        	$$h{$tt}[$j]=log2((($$h{$tt}[$j]+$pseudocount/4)/($sum+$pseudocount))/$background[$tt]);
        	#print "$$h{$tt}[$j]\t";
        }
        #print "\n";
    }
    #print "$h\n";
    return $h;
    
    
}

sub log2 {
  my $n = shift;
  return log($n)/log(2);
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
sub median{
    my ($arr)=@_;
    my @brr=sort{$a<=>$b} @$arr;
    my $n=@brr;
    if($n==0){
        return 0;
    }

    else{

        if($n%2==0){
            return $brr[$n/2]+$brr[($n/2)-1];
        }
        else{
            return $brr[($n-1)/2];
        }
    }
}


sub sum{
    my ($arr)=@_;
    my $summ=0;
    foreach my $t (@$arr){
        $summ+=$t;
    }
    return $summ;
}

sub mean{
    my ($arr)=@_;
    my $n=@$arr;
    if($n==0){
    	return 0; 
    }
    else{
	    my $summ=0;
	    foreach my $t (@$arr){
	        $summ+=$t;
	    }
	    my $n=@$arr;
	    return int($summ*1000000/$n)/100;
	  }
}

sub min{
	my ($a)=@_;
	my $min=$$a[0];
	foreach my $s (@$a){
		if($min>$s){
			$min=$s;
		}
	}
	return $min;
}

sub max{
	my ($a)=@_;
	my $max=$$a[0];
	#print "$max=1\n";
	foreach my $s (@$a){
		if($max<$s){
		    #print "max:\t$max\t$s\n";
			$max=$s;
		}
	}
	return $max;
}


sub pearson{
        my ($x,$y)=@_;
        my @ax=@$x;
        my @ay=@$y;
       
		if(@ax!=@ay){
				return "error";
		}
		else{
				my $n=@ax;
				my $x2=0;my $y2=0;my $xy=0;my $s_x=0;my $s_y=0;
				for(my $i=0;$i<$n;$i++) {
						$x2 +=$ax[$i]*$ax[$i];
						$y2 +=$ay[$i]*$ay[$i];
						$xy +=$ax[$i]*$ay[$i];
						$s_x+=$ax[$i];
						$s_y+=$ay[$i];
				}
						my $pr=0;
						$pr=($n*$xy-$s_x*$s_y)/(sqrt($n*$x2-$s_x*$s_x)*sqrt($n*$y2-$s_y*$s_y)) if(($n*$x2-$s_x*$s_x)*($n*$y2-$s_y*$s_y)>0);
						#print "$pr\n";
						return "$pr";
		}
       
}



sub id_conversion{
	my ($file,$tas)=@_;
	open(IN,"$file")||die;
	my %res=();
	while(my $t=<IN>){
		next if($t=~/^\#/);
		chomp $t;
		my @a=split(/\t/,$t);
		next if($a[0] ne $tas);
		$res{"id2sym"}{$a[1]}=$a[2];
		$res{"sym2id"}{$a[2]}=$a[1];
		next if($a[4] eq "-");
		my @aa=split(/\|/,$a[4]);
		foreach my $g (@aa){
			$res{"sym2id"}{$g}=$a[1];
		}
	}
	return \%res;
	
	
}


sub permutation {
   my (@array) = @_;
   my $length =  scalar @array;
   if($length == 1){
	   #如果数组就一个元素，返回的结果集合里面只有一个单元素数组
	   return [\@array];
   }else {
	   #取出最后一个元素
	   my $tail = pop(@array);
	   #求解N-1个数字的全排列
	   my $sub_results = permutation(@array);
	   my $new_results = [];
	   foreach my $sub_result (@$sub_results){
		  my $pos = 0;
		  #取出的最后的那个元素可以插在0,1..$length-1等处
		  while($pos < $length){
		   my @new_result = @$sub_result;
		   #插入尾元素，形成新的数组
		   splice(@new_result, $pos, 0, $tail);
		   $pos++;
		   #收集结果
		   push @$new_results, \@new_result;
		  }
	   }
	   return $new_results;
   }
}



sub write_gtf{
	my ($output,$new_trans)=@_;
	open(OUT,">$output")||die;
	foreach my $c (sort keys %$new_trans){
		foreach my $g (sort keys %{$$new_trans{$c}}){
			my $index=0;
			foreach my $trs (sort keys %{$n$ew_trans{$c}{$g}}){
				my @info=split(/\t/,$trs);
				foreach my $e_c (sort keys %{$$new_trans{$c}{$g}{$trs}}){
					my $tid=$g.".T".$index;
					#print "$trs\n$e_c\n\n" if($g eq "AN1376");
					print OUT "$c\tStringTie\ttranscript\t$info[0]\t$info[1]\t.\t$info[2]\t.\tgene_id \"$g\"; transcript_id \"$tid\";\n";
					my @aa=split(/\t/,$e_c);
					my $cc=1;
					for(my $i=0;$i<@aa;$i=$i+2){
						print OUT "$c\tStringTie\texon\t$aa[$i]\t".$aa[$i+1]."\t.\t$info[2]\t.\tgene_id \"$g\"; transcript_id \"$tid\"; exon_number \"$cc\";\n";
						$cc++;
					}
					$index++;
				}
			
			}
		}
	}
	close OUT;
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




