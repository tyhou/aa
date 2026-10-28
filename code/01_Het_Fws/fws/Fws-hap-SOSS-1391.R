#BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)
library(SeqArray)
#To convert a tabixed gzipped VCF file to the GDS format use the SeqArray package.
seqVCF2GDS("/data/01_Het_Fws/fws/F1SOSS_hap_overlap_1391.vcf",
           "/data/01_Het_Fws/fws/F1SOSS_hap_overlap_1391.gds")

#To read the GDS file into R, use the seqOpen function and to check that its output matches the orginal VCF file use seqSummary.
my_vcf <-seqOpen("/data/01_Het_Fws/fws/F1SOSS_hap_overlap_1391.gds")
seqSummary(my_vcf)

# save sample identifiers
sample.id <- seqGetData(my_vcf, "sample.id")
library(moimix)
# get genomic coordinates of all variants
coords <- getCoordinates(my_vcf)
head(coords)

#Estimating the BAF matrix
#The first step in our clonality analysis is to estiamte the B-allele frequencies for each isolate directly from the sequencing depth.
isolate_baf <- bafMatrix(my_vcf)
class(isolate_baf)
str(isolate_baf)

#Estimating MOI with Fws
#An alternative way of assessing clonality is to estimte the Fws statistic. An Fws<0.95 is indicative of a clonal infection.
fws_all <- getFws(my_vcf)
hist(fws_all)
# see if our sample that we estimated is multiclonal according to fws
fws_all["W455"] < 0.95
aa <- as.data.frame(fws_all[])
write.table(aa,file = "/results/01_Het_Fws/fws/Fws_SOSS_hap_overlap_1391.txt")
