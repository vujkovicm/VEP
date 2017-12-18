requested = read.table("../input/gene.list", F, stringsAsFactors = F)
retrieved = read.table("../output/gene.bed", F, stringsAsFactors = F)

shortlist = retrieved[which(retrieved$V4 %in% requested$V1),]

requested$X = ifelse(requested$V1 %in% shortlist$V4, 0, 1)

extra = requested[which(requested$X ==1), "V1"]

# Note: if there are duplicates in the requested, then shortlist will be unique
# e.g. rows may not always overlap but that is fine

# some type of check
is.null(dim(extra))
{
	write.table(shortlist, "../output/gene.bed", row.names = F, col.names = F, quote = F, sep = "\t")
}
