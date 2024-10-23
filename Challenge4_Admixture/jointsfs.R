
namess = c("PairwiseVCFs/1_1.recode.vcf", 
           "PairwiseVCFs/1_2.recode.vcf", 
           "PairwiseVCFs/1_3.recode.vcf", 
           "PairwiseVCFs/1_4.recode.vcf", 
           "PairwiseVCFs/1_5.recode.vcf", 
           "PairwiseVCFs/2_1.recode.vcf", 
           "PairwiseVCFs/2_2.recode.vcf", 
           "PairwiseVCFs/2_3.recode.vcf", 
           "PairwiseVCFs/2_4.recode.vcf", 
           "PairwiseVCFs/2_5.recode.vcf")
modern = c(20,20,20,20,20,16,16,16,16,16)
ancient = c(3,2,2,2,1,3,2,2,2,1)
modernname = c(1,1,1,1,1,2,2,2,2,2)
ancientname = c(1,2,3,4,5,1,2,3,4,5)

for (indexxxx in 1:10) {
  FILE = readLines(namess[indexxxx])
  
MATRIX = matrix(0, nrow = 2 * modern[indexxxx] + 1, ncol = 2 * ancient[indexxxx] + 1)

for (snpline in 10:length(FILE)) {
SNP = strsplit(FILE[snpline], "\t")[[1]]
Mainland = c()
for (i in (9+1):(9 + modern[indexxxx])) {
  Mainland = c(Mainland, as.numeric(strsplit(SNP[i], "")[[1]][1]) )
  Mainland = c(Mainland, as.numeric(strsplit(SNP[i], "")[[1]][3]) )
}
mainlandsum = sum(Mainland)
Island = c()
for (i in (9 + modern[indexxxx] + 1):(9 + modern[indexxxx]  + ancient[indexxxx])) {
  Island = c(Island, as.numeric(strsplit(SNP[i], "")[[1]][1]) )
  Island = c(Island, as.numeric(strsplit(SNP[i], "")[[1]][3]) )
}
islandsum = sum(Island)
MATRIX[mainlandsum + 1, islandsum + 1] = MATRIX[mainlandsum + 1, islandsum + 1] + 1

}
MATRIX[1,1] = 250000000 - sum(MATRIX)
write.table(MATRIX, paste0(modernname[indexxxx] , "_", ancientname[indexxxx], ".csv"  )
            , row.names = F, col.names = F, quote= F)

}
