FILE = readLines("small.recode.vcf")
MATRIX = matrix(0, nrow = 5, ncol = 5)

for (snpline in 10:length(FILE)) {
SNP = strsplit(FILE[snpline], "\t")[[1]]
Mainland = c()
for (i in (9+1):(9 + 2)) {
  Mainland = c(Mainland, as.numeric(strsplit(SNP[i], "")[[1]][1]) )
  Mainland = c(Mainland, as.numeric(strsplit(SNP[i], "")[[1]][3]) )
}
mainlandsum = sum(Mainland)
Island = c()
for (i in (9+3):(9+4)) {
  Island = c(Island, as.numeric(strsplit(SNP[i], "")[[1]][1]) )
  Island = c(Island, as.numeric(strsplit(SNP[i], "")[[1]][3]) )
}
islandsum = sum(Island)
MATRIX[mainlandsum + 1, islandsum + 1] = MATRIX[mainlandsum + 1, islandsum + 1] + 1

}
MATRIX[1,1] = 100000000 - sum(MATRIX)
write.table(MATRIX, "small_joint.csv", row.names = F, col.names = F, quote= F)
plot(colSums(MATRIX)[-1])
