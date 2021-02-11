if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    BiocManager::install("GenomicRanges")
}
library(GenomicRanges, warn.conflicts=F);
makeGRangesFromEncodeNarrowPeak <- function(filename) {
    narrowPeak_df = read.table(filename, 
           col.names=c('chrom', 
                       'comprehensive_start', 
                       'comprehensive_end',   
                       'name', 'score', 'strand',  
                       'signal',
                       'minus_log_p',     
                       'minus_log_q','offset'))
    narrowPeak_df$score = NULL; # always 0
    narrowPeak_df$strand = NULL; # always '.'
    narrowPeak_df$name = NULL; # always '.'
    narrowPeak_df$minus_log_p = NULL; # always -1
    
    narrowPeak_gr = makeGRangesFromDataFrame(narrowPeak_df, starts.in.df.are.0based=T, 
                                             keep.extra.columns=T, 
                                             start.field="comprehensive_start",
                                             end.field="comprehensive_end",
                                             ignore.strand=T);
    return(narrowPeak_gr);
}
makeGRangesFromNarrowPeak <- function(filename) {
    narrowPeak_df = read.table(filename, 
        col.names=c('chrom', 
                 'comprehensive_start', 
                 'comprehensive_end',   
                 'name', 'score', 'strand',  
                 'signal', 'skip1', 'skip2', 
                 'overall_summit_off', 
                 'minus_log_p',     
                 'minus_log_q',     
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
    narrowPeak_df$name = NULL;

    narrowPeak_gr = makeGRangesFromDataFrame(narrowPeak_df, starts.in.df.are.0based=T, 
                                                keep.extra.columns=T, 
                                                start.field="comprehensive_start",
                                                end.field="comprehensive_end",
                                                ignore.strand=T);
    return(narrowPeak_gr);

}

