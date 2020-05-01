# is the destination path there?
destPathDir = "/Volumes/onishlab_shared/PROJECTS/04_elt2_ChIPseq_Erin_David/"
NASpath=dir(destPathDir)
stopifnot(length(NASpath) > 0)

library(XLConnect)

xlsFilename="annotatedPeaks.xlsx"
annotatedPeaksWorkbook = loadWorkbook(str_c(destPath, xlsFilename, sep="/"), create=T)

# is the data file here?
filepath=dir(".", pattern="annotatedPeaks.rds")
annotatedPeaks = readRDS(filepath)
apmdata = as.data.frame(mcols(annotatedPeaks))
chroms = as.character(seqnames(annotatedPeaks))
starts = start(annotatedPeaks)
ends = end(annotatedPeaks)
rownames(df) = str_c(df$name, df$feature, sep=".")
df = data.frame(aname=rownames(df),chrom=chroms,start=starts,end=ends,apmdata)
df$fromOverlappingOrNearest = NULL
df$variance = length(apply(apmdata[,c('LE_nonNormed','L1_nonNormed','L3_nonNormed')], 1, var))
sheetname = format(Sys.time(), "%y-%m-%d-%H-%M")
createSheet(annotatedPeaksWorkbook, sheetname)
writeWorksheet(annotatedPeaksWorkbook, df, sheetname)
saveWorkbook(annotatedPeaksWorkbook)
