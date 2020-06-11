# closest_gene_calling_v10.pl
calling the closest candidate target gene for Chip-seq peaks 


Options:

       -gff                                   * gff file 
       -peak_file                             * ur input peak bed files, "*chr_name  *start  *end  length  summit  .."
       -summit_column                         Int. which column is summit info in your peak_file, count from 1. Default is unset
                                              without setting summit_column, script will use middle site(mid) of peak to call gene
                                          
       -summit_shift                          Int. when summit/mid is on gene body and very close to TSS,
                                              shift several bp to shift summit/mid out of gene to promoter, Default is 0.
                                              For TF ChIP, summit_shift is better not too long, less than 20 suggested.
                                          

       -calling_gene_by_closest               Flag. call closest one gene from peak summit or mid(cis regulated).Default is TRUE
       -abcd                                  Flag. required by -calling_gene_by_closest,call closest genes cis and trans.

       -thread_limit                          Int. thread_limit, how many threads want to be used, Default is 4
       -distance                              Int. distance, required by -calling_gene_by_distance and -calling_gene_by_both, default is 2000 (bp)



       -chromosome_name_file                  chromosome_name_file,if chr_name in gff is diff with ur chr_name in peak_file,
                                              plz set it, "chr_name(in gff)  chr_name(in peak_file)"
                                          
       -common                                Flag. if want to output gene common name
       -cds                                   Flag. if want to use cds coordinate as gene coordinate

                                              * means required.
                                          
       -help|h            brief help message
       
       
# abcd means four types of binding style(position relationship between TF and its two closest genes of each side), a, two genes are same direction of positive. c, same direction of negative. b, left side gene is positive direction and right side gene is negative direction. d, left gene is negative direction and right side is positive direction.




# zqWinSGR-v4.pl
=head1 NAME

script - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

script [options] [args ...]

Options: 

	-feature_file              * Feature_file or chromosome coordinate file, format:"chr  start  end"
	-socre_file                * socre_file, like .sgr .bed(counted)

	-chrom_column              Int. In feature_file, which column is chromosome name, counting from 1, Default is 1
	-start_colum               Int. In feature_file, which column is start site, counting from 1, Default is 2
	-end_column                Int. In feature_file, which column is end site, counting from 1, Default is 3
	-direction_column|e        Int. In feature_file, which column is direction, counting from 1, Default is 4
	
	-start_range               Defined which range be scanned, Default is -1000,0
	-end_range                 Defined which range be scanned, Default is 0,1000
	
	-start_bin                 Int. Window length for scanning start_range, units is bp, Default is 50 
	-end_bin                   Int. Window length for scanning end_range, units is bp, Default is 50   

	
	#for relative range
	-bin_count|b               Int. How many bins you want to check, Default is 20   
	-merge                     Bool. whether merge start, end, and rel score together, Default is False
	-read_count|rc             Int. read_count for normalization the score , Default is not to do normalization.
	
	#for output
	-output_folder|of          Output result folder, Default is same folder with feature_file
	-outout_name|on            Output result name, Default is same feature_file name combine socre_file name
	

	
	-help			           brief help message
