library(DEWSeq)
library(IHW)
library(tidyverse)
library(data.table)
library(ggrepel)

count_matrix <- fread('Data/Counts/Counts_matrix_w100s20.txt.gz',stringsAsFactors = F,sep = "\t", header = T)
count_matrix <- column_to_rownames(count_matrix,'unique_id')
col_data <- data.frame(type=c('Gal','Gal','Lo','Lo','Norm','Norm'),row.names = colnames(count_matrix))

annotation_file <- 'Data/Genome/SLWD_w100s20.txt.gz'

ddw <- DESeqDataSetFromSlidingWindows(countData = count_matrix,colData = col_data,annotObj = annotation_file,design = ~type)

ddw <- estimateSizeFactors(ddw)
sizeFactors(ddw)


keep <- rowSums(counts(ddw)) >= 10
ddw <- ddw[keep,]

#
rm(keep)
gc()
#


ddw_mRNAs <- ddw[rowData(ddw)[,"gene_type"]=="protein_coding",]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)
sizeFactors(ddw)



ddw_tmp <- ddw
ddw_tmp <- estimateDispersions(ddw_tmp,fitType="local",quiet=TRUE)
ddw_tmp <- nbinomWaldTest(ddw_tmp)

tmp_significant_windows <- 
  results(ddw_tmp,
          contrast = c("type", "Norm", "Lo"),
          tidy = TRUE,
          filterFun = ihw) %>% 
  dplyr::filter(padj < 0.05) %>% 
  .[["row"]]
rm("ddw_tmp")
gc()


ddw_mRNAs <- ddw_mRNAs[ !rownames(ddw_mRNAs) %in% tmp_significant_windows, ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs)
sizeFactors(ddw) <- sizeFactors(ddw_mRNAs)
rm(list = c("tmp_significant_windows","ddw_mRNAs"))
gc()
sizeFactors(ddw)


decide_fit <- TRUE

parametric_ddw <- estimateDispersions(ddw,fitType='parametric')
if(decide_fit){
  local_ddw <- estimateDispersions(ddw,fitType="local")
}

plotDispEsts(parametric_ddw,main="Parametric Fit")



if(decide_fit){
  plotDispEsts(local_ddw, main="Local fit")
}






parametricResid <- na.omit(with(mcols(parametric_ddw),abs(log(dispGeneEst)-log(dispFit))))
if(decide_fit){
  localResid <- na.omit(with(mcols(local_ddw),abs(log(dispGeneEst)-log(dispFit))))
  residDf <- data.frame(residuals=c(parametricResid,localResid),
                        fitType=c(rep("parametric",length(parametricResid)),
                                  rep("local",length(localResid))))
  summary(residDf)
}


if(decide_fit){
  ggplot(residDf, aes(x = residuals, fill = fitType)) + 
    scale_fill_manual(values = c("darkred", "darkblue")) + 
    geom_histogram(alpha = 0.5, position='identity', bins = 100) + theme_bw()
}

summary(parametricResid)
if(decide_fit){
  summary(localResid)
  if (median(localResid) <= median(parametricResid)){
    cat("chosen fitType: local")
    ddw <- local_ddw
  }else{
    cat("chosen fitType: parametric")
    ddw <- parametric_ddw
  }
  rm(local_ddw,parametric_ddw,residDf,parametricResid,localResid)
}else{
  ddw <- parametric_ddw
  rm(parametric_ddw)
}


ddw <- estimateDispersions(ddw, fitType = "local", quiet = TRUE)
ddw <- nbinomWaldTest(ddw)
plotDispEsts(ddw)


#############################33
#Get results
resultWindows <- resultsDEWSeq(ddw,
                               contrast = c("type", "Norm", "Lo"),
                               tidy = TRUE) %>% as_tibble

resultWindows

resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                                               data = resultWindows,
                                               alpha = 0.05,
                                               nfolds = 10))

result.fixed <- fdrtool::fdrtool(resultWindows$stat,statistic = "normal")
resultWindows$qvalue <- p.adjust(result.fixed$pval,method = "BH")


resultWindows <- resultWindows %>% 
  mutate(significant = resultWindows$p_adj_IHW < 0.1)

sum(resultWindows$significant)

resultWindows %>%
  filter(significant) %>% 
  arrange(desc(log2FoldChange)) %>% 
  .[["gene_name"]] %>% 
  unique %>% 
  head(40)
colnames(resultWindows)[21] <- "padj"
resultRegions <- extractRegions(windowRes  = resultWindows,
                                padjCol    = "padj",
                                padjThresh = 0.1, 
                                log2FoldChangeThresh = 0.5) %>% as_tibble

resultRegions$chromosome <- sub("^","chr",resultRegions$chromosome)

toBED(windowRes = resultWindows,
      regionRes = resultRegions,
      fileName  = "Data/Counts/Norm_vs_Gal_enrichedWindowsRegions.bed")
