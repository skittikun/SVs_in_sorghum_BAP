

dt <- "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/terra_sort_merge_CNV_352.vcf"
dat <- read.table(dt, quote = "", fill = T, header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F)#, skip = (startline(dt)-1))
dat$INFO <- gsub(";IMPRECISE", ";IMPRECISE=1", gsub(";SECONDARY", ";SECONDARY=1", dat$INFO))
















