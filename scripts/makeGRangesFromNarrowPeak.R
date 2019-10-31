if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    BiocManager::install("GenomicRanges")
}
library(GenomicRanges);
makeGRangesFromNarrowPeak <- function(filename) {
    narrowPeak_df = read.table(filename, 
        col.names=c('chrom', 
                 'comprehensive_start', 
                 'comprehensive_end',   
                 'strand', 'score', 'name',  
                 'signal', 'skip1', 'skip2', 
                 'overall_summit_off', 
                 'minus_log_q',     
                 'minus_log_p',     
                 'peak_1_start',    
                 'peak_1_end',      
                 'peak_1_signal',   
                 'peak_1_summit',   
                 'peak_2_start',    
                 'peak_2_end',      
                 'peak_2_signal',   
                 'peak_2_summit'    
                 ));

    narrowPeak_df$skip1 = NULL;
    narrowPeak_df$skip2 = NULL;
    narrowPeak_df$strand = NULL;

    narrowPeak_gr = makeGRangesFromDataFrame(narrowPeak_df, starts.in.df.are.0based=T, 
                                                keep.extra.columns=T, 
                                                start.field="comprehensive_start",
                                                end.field="comprehensive_end",
                                                ignore.strand=T);
    return(narrowPeak_gr);

}

