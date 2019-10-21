## read in real data (assume individual vcfs ##
###############################################

# Get a list of VCFs and sample names from the vcfs
vcf_list = list.files(path=vcf_dir,pattern = "\\.vcf$")
sample_ID_list = sapply(vcf_list, USE.NAMES = FALSE, function(x){strsplit(x, ".vcf", fixed=T)[[1]][1]} )

# Ths will contain a list of vcfs and samples processed.   
# It should match vcf_list and sample_ID_list, unless there is an issue
#vcf_list_new = {} # this will contain those with non-empty data
sample_ID_list_new = {}

# Read in all the VCFs to test_data.  This is basically a merged VCF
for (i in 1:length(vcf_list)) {

    message("[", i, "/", length(vcf_list), "] Reading ", paste(vcf_dir,vcf_list[i],sep=""))
    tmp = read.table(file=paste(vcf_dir,vcf_list[i],sep=""), header=T)
    
    colnames(tmp) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_ID_list[i])

    tmp$CHROM = as.character(tmp$CHROM)
    tmp$CHROM[!is.na(as.numeric(tmp$CHROM))] = paste("chr",tmp$CHROM[!is.na(as.numeric(tmp$CHROM))],sep="") # some data have CHROM as a number e.g. 16 instead of chr16

    tmp$REF = as.character(tmp$REF)
    tmp$ALT = as.character(tmp$ALT)
    
    tmp = tmp[,c("CHROM",	"POS", "REF", "ALT", "INFO", sample_ID_list[i])]
    
    tmp$INFO = sapply(tmp$INFO, function(x){strsplit(stringr::str_extract(as.character(x),"DP=[0-9]+"),"=")[[1]][2]}) # get read depth of each locus
    tmp$DP = as.numeric(tmp$INFO)
    tmp$INFO = NULL
    
    tmp = subset(tmp,nchar(tmp$REF)==1 & nchar(tmp$ALT)==1)
    tmp = subset(tmp, tmp$DP>100)
    tmp$DP = NULL
    
    if (exists("test_data")) {
      test_data = merge(test_data,tmp,by=c("CHROM","POS", "REF"),all=TRUE)
      test_data = test_data %>% mutate(ALT = coalesce(ALT.x,ALT.y)) %>% select(-c("ALT.x","ALT.y")) # assume only one possible valid ALT
      test_data = test_data[,c("CHROM","POS","REF","ALT",colnames(test_data)[which(!colnames(test_data)%in%c("CHROM","POS","REF","ALT"))])]
    } else {
      test_data = tmp
    }
    
    rm(tmp)
    
    #vcf_list_new = c(vcf_list_new, vcf_list[i])
    sample_ID_list_new = c(sample_ID_list_new, sample_ID_list[i])
}
rm(i)
rm(vcf_list)
rm(sample_ID_list)

# Read in sampleToSubject map
sampleTosubject =  read.table(file=paste("sampleToSubject.txt",sep=""), header=T, stringsAsFactors =F) # nrow could be > number of samples in test_data
sampleTosubject = subset(sampleTosubject, sampleTosubject$SAMPLE%in%sample_ID_list_new)
sampleTosubject = sampleTosubject[match(sample_ID_list_new, sampleTosubject$SAMPLE),]
colnames(sampleTosubject) = c("sample_ID","subject_ID")
sampleTosubject = sampleTosubject[order(sampleTosubject$subject_ID),]

for (i in seq(length(sampleTosubject$sample_ID))) {
  message("[", i, "/", length(sampleTosubject$sample_ID), "] Converting ", sampleTosubject$sample_ID[i])
  
  # Save only the genotype (0/0, 0/1, 1/1)
  test_data[,sampleTosubject$sample_ID[i]] = sapply(test_data[,sampleTosubject$sample_ID[i]], function(x){strsplit(as.character(x),":")[[1]][1]})
  # Convert genotype to base (A/C, etc) if there is a genotype, else NA/NA
  test_data[,sampleTosubject$sample_ID[i]]  = mapply(
    FUN=function(x,ref,alt){
      genotype_num = strsplit(as.character(x),"/")[[1]]
      genotype_letter = c("NA","NA")
      genotype_letter[1] = ifelse(genotype_num[1]==".","NA",  ifelse(genotype_num[1]=="0",ref,alt))
      genotype_letter[2] = ifelse(genotype_num[2]==".","NA",  ifelse(genotype_num[2]=="0",ref,alt))
      paste(genotype_letter[1],"/",genotype_letter[2],sep="")
    }
    ,test_data[,sampleTosubject$sample_ID[i]], test_data$REF,test_data$ALT
  )
  # If genotype is missing (NA/NA), just make it NA
  test_data[,sampleTosubject$sample_ID[i]][which(test_data[,sampleTosubject$sample_ID[i]]=="NA/NA")] = NA
}
rm(i)

# Select only for loci with alt allele
test_data = subset(test_data,test_data$ALT!=".")

## subset test data ##

tmp=apply(test_data[,5:ncol(test_data)], 2, function(x){sum(table(x))}) # number loci recorded for each sample
sample_ID_tmp = names(tmp[which(tmp>20)])

test_data.subset = test_data[,c(colnames(test_data)[1:4],sample_ID_tmp)]
sampleTosubject.subset = subset(sampleTosubject,sampleTosubject$sample_ID%in%sample_ID_tmp)
