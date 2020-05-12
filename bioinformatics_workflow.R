#Reading GPL Files downloaded from NCBI GEO

anno <- function(file){
  
  #Read in the gpl files that have created previously using Biomart, 
  #Each dataset contain the probe, HGNC and Entrez ID  
  
  df1 <- read.csv(file)
  df1 <- as.data.frame(df1)
  df1 <- df1[, -1]
  
  colnames(df1) <- c('ProbeID', 'Entrez', 'Symbol')
  
  return(df1)
  
}


#Reading .CEL files, merging and normalising 

preprocess <- function(files, cols, name){
  
  library(oligo)
  library(oligoClasses)
  library(BiocParallel)
  library(preprocessCore)
  library(robustHD)
  library(lumi)
  library(limma)
  library(gcrma)
  
  samples <- list()
  colno <- list()
  cc <- c()
  
  for (file in files){
    
    read <- list.celfiles(paste('file_directory', 
                                file, sep=''), full.name=TRUE)
    read <- read.celfiles(read)
    read <- oligo::rma(read)
    read <- exprs(read)
    
    cc <- append(cc, colnames(read))
    colno <- append(colno, ncol(read))
    samples <- append(samples, list(read))
    
  }
  
  rows <- rownames(samples[[1]])
  dfnew <- data.frame(row.names=rows)
  
  for (idx in samples){
    
    dfnew <- cbind(dfnew, idx)
    
  }
  
  logX <- as.matrix(dfnew)
  exprs <- normalize.quantiles(logX, copy=TRUE)
  expr <- as.data.frame(exprs)

  rr <- rownames(dfnew)
  cc <- colnames(dfnew)
  
  colnames(expr) <- cols
  rownames(expr) <- rr
  
  return(list(expr, logX))
  
  
}


#Differential Gene expression using Limma

sigde <- function(input, name, anno, fc, colna){
    
  library(limma)
  library(data.table)
  library(vsn)
  library(dplyr)
  
  design <- model.matrix(~0 + factor(c(rep(1, 3), rep(2, 3), rep(3, 3))))
  colnames(design) <- c("Case", "Control", "Normal")
  
  fit <- lmFit(input, design, method="ls")
  contrast.matrix <- makeContrasts(Case-Normal, 
                                   Case-Normal,
                                   Normal-Control,
                                   levels=design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
  
  results <- decideTests(fit2, method="separate", 
                         adjust.method="BH", 
                         p.value=0.05, lfc=2)
  
    
  dfList_up <- list()
  dfList_dwn <- list()
  
sigs <- list()
  
  for (k in 1:2){
    
    sig <- topTable(fit2,
                    n=Inf,
                    adjust="BH",
                    coef=k,
                    sort.by="P",
                    p.value=0.05,
                    lfc=fc)
    
 
    up <- sig[which(sig$logFC > 1), ]
    dwn <- sig[which(sig$logFC < 1), ]
    
    dfList_up <- append(dfList_up, list(up))
    dfList_dwn <- append(dfList_dwn, list(dwn))
    sigs <- append(sigs, list(sig))
    
  }
  
  sig_genes <- list(dfList_up, dfList_dwn)
  genes <- list()

  for (r in 1:2){
    
    int_genes <- Reduce(intersect, lapply(sig_genes[[r]], rownames))
    df_data <- input[int_genes, ]
    df_data <- na.omit(df_data)
    genes <- append(genes, list(df_data))
    
  }
  
  gene_names <- list('Up-reguated', 'Down-regulated')
  output2 <- rbind(genes[[1]], genes[[2]])
  output <- genes[[1]]
  
  all <- list(output2, output)
  out <- list()
  
  for (p in 1:2){
    
    output <- data.frame(all[[p]])
    anno <- data.frame(anno)
    
    x_input <- setDT(output, keep.rownames = TRUE)[]
    x_input <- data.frame(x_input)
    
    colnames(output)[1] <- c("probeID")
    colnames(anno)[1] <- c("probeID")
    
    dataset_n <- merge(anno, output, by='probeID')
    dataset_n <- dataset_n[, -1]
    
    dataset_n <- dataset_n[!duplicated(dataset_n$Symbol), ]
    dataset_n <- na.omit(dataset_n)
    rownames(dataset_n) <- dataset_n$Symbol
    
    dataset_n <- as.data.frame(dataset_n)
    dataset_n <- select(dataset_n,-c(Entrez,Symbol))
    
    dataset_n <- as.matrix(dataset_n)
    colnames(dataset_n) <- colna
    rownames(dataset_n) <- rownames(dataset_n)
    
    out <- append(out, list(dataset_n))
  }
    
  chge <- list('up','down')
  
  for (j in 1:2){
    
    logX <- data.frame(genes[[j]])
    anno <- data.frame(anno)
    
    x_input <- setDT(logX, keep.rownames = TRUE)[]
    x_input <- data.frame(x_input)
    
    colnames(x_input)[1] <- c("probeID")
    colnames(anno)[1] <- c("probeID")
    
    dataset_new <- merge(anno, x_input, by='probeID')
    dataset_new <- dataset_new[, -1]
    
    dataset_new <- dataset_new[!duplicated(dataset_new$Symbol), ]
    dataset_new <- na.omit(dataset_new)
    rownames(dataset_new) <- dataset_new$Symbol
    
    dataset_new <- as.data.frame(dataset_new)
    dataset_new <- select(dataset_new,-c(Entrez,Symbol))
    colnames(dataset_new) <- colna
    
  }
  
  return(out[[2]]) #1 all, 2 upregulated 
  
}



#Soft-thresholding power determination 

pwr <- function(dataset, name){
  
  library(WGCNA)
  enableWGCNAThreads(nThreads = 2)
  allowWGCNAThreads(nThreads = 4)
  options(stringsAsFactors = FALSE) 
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(dataset, powerVector = powers, verbose = 5)

  
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
       labels=powers,cex=0.9, col="red")
  abline(h=0.85,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
       labels=powers, cex=0.9, col="red")
  
}

#Constructing co-expression network

TOM <- function(dataset, pwr, name){
  
  library(WGCNA)
  
  enableWGCNAThreads(nThreads <- 2)
  allowWGCNAThreads(nThreads = 4)
  options(stringsAsFactors = FALSE)
  
  adj <- abs(cor(dataset, use="p"))^pwr
  dissTOM <- TOMdist(adj)
  simTOM <- TOMsimilarity(adj)
  
  #hierTOMa <- hclust(as.dist(dissTOM), method="average")
  #hierTOMc <- hclust(as.dist(dissTOM), method="complete")

  #hieradja <- hclust(as.dist(adj), method="average")
  #hieradjc <- hclust(as.dist(adj), method="complete")
 
  #hierTOMSa <- hclust(as.dist(simTOM), method="average")
  hierTOMSc <- hclust(as.dist(simTOM), method="complete")

  #HIERTOM
  PAM <- as.character(cutreeStaticColor(hierTOMSc, cutHeight=0.1, minSize=50)) 
  Gene_Modules <- labels2colors(cutreeDynamic(hierTOMSc, method="tree", cutHeight=0.98))
  HierClust_PAM <- labels2colors(cutreeDynamic(hierTOMSc, distM= simTOM , cutHeight = 0.99,
                                               deepSplit=3, pamRespectsDendro = FALSE))
  
  par(mfrow = c(1,1))
  plotDendroAndColors(hierTOMSc,
                      colors = data.frame(PAM, HierClust_PAM),
                      dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                      cex.axis = 1.2)
  
  return(list(adj, dissTOM, Gene_Modules))
  
}


#Constructing of gene eigengene 

eigenetic_network <- function(dataset, colorh1, colna, name){
  
  library(purrr)
  
  change <- t(dataset)
  colnames(change) <- colna
  
  datME <- moduleEigengenes(t(change),colorh1)$eigengenes
  MET <- orderMEs(cbind(datME))
  
  par(mfrow=c(1, 1))
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                        marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                        xLabelsAngle = 90)
    
  datKME <- signedKME(dataset, datME, outputColumnName="") 
  
  return(datKME)
  
}
 

#Module membership determination 

intramodcon <- function(dataset, adj, colorlevels, colorh1, name, km, lbs){
  
  which.color=colorlevels;
  restrictGenes=colorh1==which.color 
  verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                     (km[restrictGenes, paste(which.color, sep="")])^3, 
                     col=which.color,
                     xlab="Intramodular Connectivity",
                     ylab="Module Membership")

}


#Functional enrichment analysis for each gene module 

enrichment <- function(dataset, colorh1, datKME, anno, name){
  
  library(rowr)
  library(dplyr)

  y <- y[, -1]

  intModules <- table(colorh1)
  intModules <- as.data.frame(intModules)
  intModules <-intModules$colorh1
  intModules <- as.character(intModules)
  
  dat <- data.frame() 
  colrs <- c()
  
  for (color in intModules){
    
    FilterGenes <- abs(subset(datKME, select=c(color))) > 0.40
    
    genes <- dimnames(data.frame(dataset))[[2]][FilterGenes]
    entr <- dimnames(data.frame(dataset))[[2]][FilterGenes]
    
    entr <- filter(anno, Symbol %in% entr)
    dat <- rowr::cbind.fill(dat, entr$Entrez, fill = NA)
    colrs <- append(color, colrs)

    write.csv(genes, file=paste(color, '.csv', sep=''))

    if (is_empty(genes) != TRUE){
      
      library(enrichR)
      
      dbs <- c("GO_Biological_Process_2018", "KEGG_2019_Human", 
               "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
      enriched <- enrichr(as.character(genes), dbs)
  
      
      EnrichRBP <- enriched[['GO_Biological_Process_2018']]
      EnrichRMF <- enriched[['GO_Molecular_Function_2018']]
      EnrichRCC <- enriched[['GO_Cellular_Component_2018']]
      EnrichRKE <- enriched[['KEGG_2019_Human']]
      
      enriched_datasets <- list(EnrichRBP, EnrichRMF, EnrichRCC, EnrichRKE)
      enriched_labels <- c('GO_Biological_Process_2018', 'GO_Molecular_Function_2018',
                           'GO_Cellular_Component_2018', 'KEGG_2019_Human')
      
      for (read in 1:4){
        
        enrichset <- enriched_datasets[[read]][which(enriched_datasets[[read]]$Adjusted.P.value < 0.9), ]
        
        if (is_empty(enrichset) != TRUE & (nrow(enrichset) > 0)){
          write.csv(enrichset, file = paste(color, enriched_labels[[read]], ".csv", sep = ""))
          
        }
      }
    }
  }
  
  dat <- dat[,-1]
  colnames(dat) <- colrs
  return(dat)
  
}

#Hub gene determination 

hub_genes <- function(dataset, colorh1, pwr){
  
  hub_genes <- chooseTopHubInEachModule(dataset, colorh1, 
                                        omitColors = "grey", 
                                        power = pwr, 
                                        type = "signed")
  
  print(as.data.frame(hub_genes))
  
}


#extracting highly connected genes from specific gene modules 

highly_connected_genes <- function(dataset, datKME, color){
  
  library(clusterProfiler)
  
  FilterGenes <- abs(subset(datKME, select=c(color))) > 0.85
  genes <- dimnames(data.frame(dataset))[[2]][FilterGenes]
  genes <- na.omit(genes)
  
  eg <- bitr(as.character(genes), fromType="SYMBOL", toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
  
}

