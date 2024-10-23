FILE = readLines("GHIST-secondary-contact.vcf")
MATRIX = matrix(0, nrow = 45, ncol = 17)

for (snpline in 10:length(FILE)) {
SNP = strsplit(FILE[snpline], "\t")[[1]]
Mainland = c()
for (i in 10:(9+22)) {
  Mainland = c(Mainland, as.numeric(strsplit(SNP[i], "")[[1]][1]) )
  Mainland = c(Mainland, as.numeric(strsplit(SNP[i], "")[[1]][3]) )
}
mainlandsum = sum(Mainland)
Island = c()
for (i in (10+22):39) {
  Island = c(Island, as.numeric(strsplit(SNP[i], "")[[1]][1]) )
  Island = c(Island, as.numeric(strsplit(SNP[i], "")[[1]][3]) )
}
islandsum = sum(Island)
MATRIX[mainlandsum + 1, islandsum + 1] = MATRIX[mainlandsum + 1, islandsum + 1] + 1

}
MATRIX[1,1] = 100000000 - sum(MATRIX)
write.table(MATRIX, "joint.csv", row.names = F, col.names = F, quote= F)
MATRIX[45,17] = .1
heatmap(log(MATRIX),Colv = NA, Rowv = NA, scale='none')
