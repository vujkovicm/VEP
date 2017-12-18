library("reshape2")

# import annotation
anno = read.table("../output/region.anno.tab", T, sep = "\t", stringsAsFactors=F)

# import counts
mi = read.table("../output/counts_MI.hwe", T, stringsAsFactors = F)
dm = read.table("../output/counts_DM.hwe", T, stringsAsFactors = F)
ldl = read.table("../output/means_LDL.qassoc.means", T, stringsAsFactors = F)
hdl = read.table("../output/means_HDL.qassoc.means", T, stringsAsFactors = F)
tg = read.table("../output/means_TG.qassoc.means", T, stringsAsFactors = F)
tc = read.table("../output/means_TC.qassoc.means", T, stringsAsFactors = F)

ldl = subset(ldl, ldl$VALUE == "MEAN")
hdl = subset(hdl, hdl$VALUE == "MEAN")
tg  = subset(tg, tg$VALUE == "MEAN")
tc  = subset(tc, tc$VALUE == "MEAN")

colnames(ldl)[4:6] = c("HOM.LDL", "HET.LDL", "WT.LDL")
colnames(hdl)[4:6] = c("HOM.HDL", "HET.HDL", "WT.HDL")
colnames(tg)[4:6] = c("HOM.TG", "HET.TG", "WT.TG")
colnames(tc)[4:6] = c("HOM.TC", "HET.TC", "WT.TC")


mi = mi[which(mi$TEST %in% c("ALL", "AFF")),]
mi[mi$TEST == "AFF", "TEST"] = "MI"
dm = dm[which(dm$TEST %in% c("AFF")),]
dm$TEST = "DM"

all = rbind(mi, dm)
all = all[,c("SNP", "TEST", "GENO")]
all.wide = reshape(all, idvar = "SNP", timevar = "TEST", direction = "wide")
all.p = mi[which(mi$TEST == "ALL"), c("SNP","P")]
colnames(all.p)[2] = "HWE.ALL"
all.wide = merge(all.wide, all.p, by = "SNP")

all.wide = merge(all.wide, tg[,c(2,4:6)], by = "SNP")
all.wide = merge(all.wide, tc[,c(2,4:6)], by = "SNP")
all.wide = merge(all.wide, ldl[,c(2,4:6)], by = "SNP")
all.wide = merge(all.wide, hdl[,c(2,4:6)], by = "SNP")

final = merge(anno, all.wide, by.x = "Variant", by.y = "SNP")

# fill in gene names
for (i in 2:(nrow(final)-1))
{
        if (final$Gene[i] == "." & final$Gene[i-1] == final$Gene[i+1])
        {
                final$Gene[i] = final$Gene[i-1]
        }
        else if (final$Gene[i] == "." & final$Position[i] == final$Position[i+1])
        {
                final$Gene[i] = final$Gene[i+1]
        }
        else if (final$Gene[i] == "." & final$Position[i] == final$Position[i-1])
        {
                final$Gene[i] = final$Gene[i-1]
        }
}

write.table(final, "../output/project_summary.tab", sep = "\t", row.names=F, col.names=T, quote = F)

