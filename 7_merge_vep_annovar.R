### the user-defined R functions
source(paste(getwd(), "/rfunc_select_the_most_severe_consequence.R", sep=""))
###=========================================================

require(data.table) ## 1.9.2

gene = read.table("../input/gene.list", F, stringsAsFactors = F)
gene = unlist(gene)

lookup = read.table("../input/consequence_severity.txt", T, stringsAsFactors = F)

a = read.table(file="../output/region.lof.out", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
a = subset(a, a$SYMBOL %in% gene)

b = read.table("../output/region.vcf.anno.avout.hg19_multianno.txt", F, sep = "\t", stringsAsFactors = F, skip = 1, quote = NULL)
c = read.table("../output/region.vcf.anno.avout.hg19_multianno.txt", F, sep = "\t", stringsAsFactors = F, nrow = 1)
c = unlist(c)
names(b)[1:82] = c
b = b[,1:82]
b = unique(b)
b$Uploaded_variation = paste(b$Chr, b$Start,b$Ref,b$Alt, sep = ":")

a = unique(a)
a$id = row.names(a)

### select the most severe consequence
a$mostSevConsq <- sapply(X = a$Consequence, FUN = selMostSev)
a = merge(a, lookup, by.x = "mostSevConsq", by.y = "LOF_Consequence", all.x = T)

sev = a[,c("id","Uploaded_variation", "Allele", "SYMBOL", "FLAGS", "TREMBL", "mostSevConsq", "Severity")]
sev$second_filter = ifelse(sev$FLAGS == "." & sev$TREMBL == ".", 0, ifelse(sev$FLAGS == "." & sev$TREMBL != ".", 1, ifelse(sev$FLAGS != "." & sev$TREMBL == ".", 2, 3)))
require(data.table) ## 1.9.2
group <- as.data.table(sev)
x = group[group[, .I[Severity == min(Severity)], by = Uploaded_variation]$V1]
y = x[x[, .I[second_filter == max(second_filter)], by = Uploaded_variation]$V1]
z = y[y[, .I[id == min(id)], by = Uploaded_variation]$V1]
z$KEEP = 1
z = as.data.frame(z)
final = merge(a, z[,c("id", "KEEP")], by = "id")
final = final[order(as.numeric(final$id)), ]

mrg = merge(final, b, by = "Uploaded_variation", all = T)

fam = read.table("../output/region.bim", F, stringsAsFactors = F, sep = "\t")
fam$tmp = as.numeric(row.names(fam))
#fam$sort = row.names(fam)
fam = merge(fam, mrg, by.x = "V2", by.y = "Uploaded_variation", all.x = T)
fam = fam[,c("V1", "V2", "V4", "V5", "V6", "mostSevConsq", "SYMBOL", "Existing_variation", "Protein_position","Amino_acids", "Codons","CLINSIG","CLNDBN", "GTEx_V6_gene", "GTEx_V6_tissue", "LoF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "SIFT", "SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "tmp")]

fam[is.na(fam)] = "."
#fam$Polyphen2_VEP = ifelse(fam$PolyPhen != ".", sapply(strsplit(fam$PolyPhen, "\\("), "[", 1), ".")
#fam$PolyPhen = NULL
#fam$Polyphen2_VEP = ifelse(fam$Polyphen2_VEP == "benign", "B", ifelse(fam$Polyphen2_VEP == "probably_damaging", "D", ifelse(fam$Polyphen2_VEP == "possibly_damaging", "P", ".")))
fam$SIFT_VEP = ifelse(fam$SIFT != ".",  sapply(strsplit(fam$SIFT, "\\("), "[", 1), ".")
fam$SIFT = NULL
fam$SIFT_VEP = ifelse(grepl("_", fam$SIFT_VEP), sapply(strsplit(fam$SIFT_VEP, "\\_"), "[", 1), ".")
fam$SIFT_VEP = ifelse(fam$SIFT_VEP == "deleterious", "D", ifelse(fam$SIFT_VEP == "tolerated", "T", "."))
fam$SIFT_pred = ifelse(fam$SIFT_VEP == "D" | fam$SIFT_pred == "D", "D", ifelse(fam$SIFT_VEP == "." & fam$SIFT_pred == ".", ".", "B"))
fam$SIFT_VEP = NULL
fam$CLINVAR = ifelse(grepl("enign", fam$CLINSIG), "B", ifelse(grepl("athogenic", fam$CLINSIG), "P",  ifelse(grepl("ncertain", fam$CLINSIG), "U", ".")))
fam$CLINPHE = ifelse(grepl("not_provided", fam$CLNDBN) | grepl("not_specified", fam$CLNDBN) , "not_specified", fam$CLNDBN)
fam$CLINSIG = fam$CLNDBN = NULL
fam$LoF = ifelse(fam$LoF == ".", ".", "LoF") 

fam$MutA = ifelse(fam$MutationAssessor_pred %in% c("H", "M"), 1, ifelse(fam$MutationAssessor_pred %in% c("L"), .5, ifelse(fam$MutationAssessor_pred %in% c("N"), 0.0, 0)))
fam$MutT = ifelse(fam$MutationTaster_pred %in% c("A", "D"), 1, ifelse(fam$MutationTaster_pred %in% c("P"), .5, ifelse(fam$MutationTaster_pred %in% c("N"), 0.0, 0)))
fam$PP2_HVAR = ifelse(fam$Polyphen2_HVAR_pred %in% c("D"), 1, ifelse(fam$Polyphen2_HVAR_pred %in% c("P"), .5, ifelse(fam$Polyphen2_HVAR_pred %in% c("B"), 0.0, 0)))
fam$PP2_HDIV = ifelse(fam$Polyphen2_HDIV_pred %in% c("D"), 1, ifelse(fam$Polyphen2_HDIV_pred %in% c("P"), .5, ifelse(fam$Polyphen2_HDIV_pred %in% c("B"), 0.0, 0)))
fam$LRT = ifelse(fam$LRT_pred == "D", 1, ifelse(fam$LRT_pred == "N", 0.0, 0))
fam$FATHMM = ifelse(fam$FATHMM_pred == "D", 1, ifelse(fam$FATHMM_pred == "T", 0.0, 0))
fam$SIFT = ifelse(fam$SIFT_pred == "D", 1, ifelse(fam$SIFT_pred == "B", 0.0, 0))
fam$CLINSIG = ifelse(fam$CLINVAR %in% c("P", "U"), 1, ifelse(fam$CLINVAR == "B", 0.5 , 0))
fam$Damage.score = fam$MutA + fam$CLINSIG + fam$PP2_HDIV + fam$LRT + fam$FATHMM

colnames(fam)[names(fam)=="GTEx_V6_gene"] <- "GTEx_Gene"
colnames(fam)[names(fam)=="GTEx_V6_tissue"] <- "GTEx_Tissue"
colnames(fam)[names(fam)=="Existing_variation"] <- "rsid"
colnames(fam)[names(fam)=="mostSevConsq"] <- "Consequence"
colnames(fam)[names(fam)=="SYMBOL"] <- "Gene"
colnames(fam)[names(fam)=="V1"] <- "Chr"
colnames(fam)[names(fam)=="V2"] <- "Variant"
colnames(fam)[names(fam)=="V4"] <- "Position"
colnames(fam)[names(fam)=="V5"] <- "Ref"
colnames(fam)[names(fam)=="V6"] <- "Alt"
colnames(fam)[names(fam)=="gnomAD_AFR_AF"] <- "MAF_AFR"
colnames(fam)[names(fam)=="gnomAD_AMR_AF"] <- "MAF_AMR"
colnames(fam)[names(fam)=="gnomAD_NFE_AF"] <- "MAF_EUR"
colnames(fam)[names(fam)=="gnomAD_EAS_AF"] <- "MAF_EAS"
colnames(fam)[names(fam)=="gnomAD_SAS_AF"] <- "MAF_SAS"

fam = fam[order(as.numeric(fam$tmp)),]
fam$tmp = NULL
fam = fam[,c(1:28)]

### output
write.table(fam, file="../output/region.anno.tab", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

