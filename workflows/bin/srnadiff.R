#!/usr/bin/env Rscript
####### Loading library

library(optparse)
library(srnadiff)
library(GenomicRanges)

####### Parse option command line
option_list <- list(
	make_option(
		c("-s", "--sampleSheet"),
		type = "character",
		default = NULL,
		help = "TSV or CSV-format sample sheet file."
	),
	make_option(
		c("-a", "--annotationFile"),
		type = "character",
		default = NULL,
		help = "GTF-format annotation file."
	),
	make_option(
		c("-f", "--feature"),
		type = "character",
		default = NULL,
		help = "feature of annotation file."
	),
	make_option(
		c("-o", "--source"),
		type = "character",
		default = NULL,
		help = "source of annotation file."
	),
	make_option(
		c("-m", "--diffMethod"),
		type = "character",
		default = "DESeq2", #edgeR #baySeq
		help = "Name of the ethod used to compute the p-values"
	), 
	make_option(
		c("-n", "--normFactors"),
		type = "character",
		default = NULL,	
		help = "One size factor for each sample in the data."
	),
	make_option(
		c("-p", "--pvalue"),
		type = "integer",
		default = 1,
		help = "Adjusted p_value of differential expressed regions."
	)
)

 
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);



####### Data preparation sample sheet, bam file and annotation file 


# Sample sheet
sample_sheet <- read.delim(file = opt$sampleSheet, header = T, sep = ",")
FileName = sample_sheet[,1]


# BAM files
#bamFiles <- opt$bamPath

# Annotation file
#annotReg <- readAnnotation(fileName = opt$annotationFile, feature = "opt$feature", source = "opt$source")

# Preparation of srnadiff object

if ( is.null(opt$annotationFile) == T ) {
	srnaExp_object <- 
		srnadiffExp(
		FileName, 
		sample_sheet,
		diffMethod = opt$diffMethod,
		normFactors = opt$normFactors) 
		
} else { 

annotReg <- readAnnotation(fileName = opt$annotationFile)

srnaExp_object <- 
		srnadiffExp(
		FileName, 
		sample_sheet, 
		annotReg,
		diffMethod = opt$diffMethod,
		normFactors = opt$normFactors) 
}

# Detecting DERs and quantifying differential expression

srnaExp <- srnadiff(
	srnaExp_object)




#Visualization of the results

gr <- regions(srnaExp, pvalue = opt$pvalue)

df <- data.frame(chr = seqnames(gr),
		starts = start(gr),
		ends = end(gr),
		names = names(gr),
		scores = -log10(mcols(gr)$padj),
		strands = strand(gr))
		
write.table(df, file = "DE_regions.bed", quote = F, sep = "\t",col.names = F, row.names = F)







	
	
	
	
	
	
	
	
	
	
	
	
