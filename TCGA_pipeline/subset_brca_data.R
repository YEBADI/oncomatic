dataSubt <- TCGAquery_subtype(tumor = BRCA)
lumA <- dataSubt[which(dataSubt$PAM50.mRNA == "Luminal A"),1]
lumB <- dataSubt[which(dataSubt$PAM50.mRNA == "Luminal B"),1]
her2 <- dataSubt[which(dataSubt$PAM50.mRNA == "HER2-enriched"),1]
basl <- dataSubt[which(dataSubt$PAM50.mRNA == "Basal-like"),1]
norml <- dataSubt[which(dataSubt$PAM50.mRNA == "Normal-like"),1]
