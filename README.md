# closest_gene_calling
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
