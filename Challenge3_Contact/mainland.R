FILE = readLines("mainland.recode.vcf")
MATRIX = matrix(-1, nrow = length(FILE) - 9, ncol = 44)
counter = 1
SNPNAMES = rep(-1, length(FILE) - 9)
for (snpline in 10:length(FILE)) {
SNP = strsplit(FILE[snpline], "\t")[[1]]
SNPNAMES[counter] = SNP[2]
minicounter = 1
for (i in 10:length(SNP)) {
  MATRIX[counter, minicounter] = as.numeric(strsplit(SNP[i], "")[[1]][1])
  minicounter = minicounter + 1
  MATRIX[counter, minicounter] = as.numeric(strsplit(SNP[i], "")[[1]][3])
  minicounter = minicounter + 1
}
counter = counter + 1
}
rownames(MATRIX) = SNPNAMES
SFS = table(rowSums(MATRIX))
write.table(SFS, "mainland.csv", row.names = F, col.names = F, quote= F)
plot(SFS)
