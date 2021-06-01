getwd()
library(tidyverse)
library(Seurat)

# 1. filter FB, DC cells : filtering 된 s1-s17 count data를 얻게 된다

## 아래와 같은 command로 17 sample 모두 읽기
# s1<- as.data.frame(read.table("GSM4430459_MS.sample1.clean.data.txt",header = T,sep = ",", row.names = 1))

## 아래와 같은 command로 전부 seurat object 생성 
# s1_seurat <- CreateSeuratObject(counts = s1 , project = "project", min.cells = 3, min.features = 200)

## 같은 status 끼리 merge (LS : 1,2,5,7 / NL : 3,11,14,15,16 / H : 4,6,8,9,10,12,13,17 ) : ls만 code를 적었으나 nl, h도 알맞은 sample들끼리 merge 해주었다.
# ls <- merge(s1_seurat, y = c(s2_seurat,s5_seurat,s7_seurat), project = "project")

#data integration & clustering analysis : nl, h 모두 동일한 방법으로 integration + clusering analysis 진행
# testlist <- SplitObject(ls, split.by = "orig.ident")
# testlist <- lapply(X = testlist, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# features <- SelectIntegrationFeatures(object.list = testlist)
# immune.anchors <-FindIntegrationAnchors(object.list = testlist, anchor.features = features)
# lsfibs <- IntegrateData(anchorset = immune.anchors)
# DefaultAssay(lsfibs) <- "integrated"
 
# lsfibs <- ScaleData(lsfibs, verbose = FALSE)
# lsfibs <- RunPCA(lsfibs, npcs = 30, verbose = FALSE)
# lsfibs <- RunTSNE(lsfibs, reduction = "pca", dims = 1:30,check_duplicates = FALSE)
# lsfibs <- FindNeighbors(lsfibs, reduction = "pca", dims = 1:30)
# lsfibs <- FindClusters(lsfibs, resolution = 0.5)
# DimPlot(lsfibs, reduction = "tsne", group.by = "orig.ident")

## marker gene으로 filtering ( FB : COL1A1, DCN / DC: CD1C) : subset에서 marker gene의 발현이 뚜렷한 cluster의 number를 적어준다. FB 뿐만 아니라 DC도 똑같은 Code로 filtering을 진행하였다.
# FeaturePlot(lsfibs, features = c("COL1A1","DCN))
# markers <- FindAllMarkers(lsfibs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# lsfiltered <- subset(lsfibs, subset=seurat_clusters %in% c(cluster numbers))

## ls, nl, h에 대해 위의 과정을 다 진행한 후 셋을 merge
# fibs <- merge(lsfiltered, y = c(nlfiltered,hfiltered) , project = "project")

## cell id 얻기 
# cells <- fibs@assays[["integrated"]]@data@Dimnames[[2]]
# write.table(cells,file = 'cells.txt',sep = '\t', row.names = T, col.names = T, quote = F)

## raw data에서 해당 cell에 관한 정보만 얻어내기 : 아래의 과정을 s1부터 s17까지 다 해준다.
# s1selc <- read.table("../cells.txt")
# id <- s1selc$V1
# filtered <- filter(s1, V1 %in% id)
# colnames(filtered)[1] <- ""
# write.table(filtered, file=gzfile("s1_fil.txt.gz"), row.names = FALSE, sep = ",", quote = FALSE)

###################################################################################

#2. Clusering analysis with filtered FB & DC cells : 아래의 모든 과정을 DC에서도 동일하게 진행해주어야 합니다.

## read filtered count data & make seurat object (from 1 through 17) ## 
s1<- as.data.frame(t(read.table("s1_fil.txt.gz",header = T,sep = ",", row.names = 1)))
s1_seurat <- CreateSeuratObject(counts = s9, project = "project", min.cells = 3, min.features = 200)

## merge all seurat objects ##
allfibs <- merge(s1_seurat, y = c(s2_seurat,s3_seurat,s4_seurat,s5_seurat,s6_seurat,s7_seurat,s8_seurat,s9_seurat,s10_seurat,s11_seurat,s12_seurat,s13_seurat,s14_seurat,s15_seurat,s16_seurat,s17_seurat) , project = "project")

## do clustering anaylsis ##
allfibs <- NormalizeData(allfibs)
allfibs <- FindVariableFeatures(allfibs, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(allfibs)
allfibs <- ScaleData(allfibs, features = all.genes)
allfibs <- RunPCA(allfibs, features = VariableFeatures(object = allfibs))
allfibs <- FindNeighbors(allfibs, dims = 1:10)
allfibs <- FindClusters(allfibs, resolution = 0.5)
allfibs <- RunTSNE(allfibs, dims = 1:10,check_duplicates = FALSE)

DimPlot(allfibs, reduction = "tsne", label = TRUE, repel = TRUE)

## plot by status ##
#sample마다 status를 cell의 수만큼 반복한 character를 만든 후 plotting
table(allfibs@meta.data$orig.ident)
sample <- c(rep("ls",1201), rep("ls",444),rep("nl",801),rep("h",2156),rep("ls",782),rep("h",922),rep("ls",147), rep("h",232),rep("h",121),rep("nl",473),rep("h",225),rep("h",62),rep("nl",188), rep("nl",10), rep("h",291))
test@meta.data$orig.ident <- sample

DimPlot(test, reduction = "tsne", group.by = "orig.ident")

###################################################################################

# 3. pyscenic : jupyter notebook 결과 파일로 첨부

## pyscenic input file ##
write.table(as.matrix(GetAssayData(object = allfibs, slot = "counts")),"counts.txt", sep = ",", quote = F, row.names = T)


# 4. cellphonedb 

## filter FB sub clusters & get count dat - COL6A5+COL18A1+ / COL11A1+LAMC3+ / APOE+ABCA+ 이 세 sub cluster에 대해 모두 동일한 code 적용
### 이 과정은 DC sub cluster filtering 및 count dat 얻는 데에 동일하게 적용해주어야 합니다. - irf8 / lamp3 / mmp12
FeaturePlot(allfibs, features = c("COL6A5","COL18A1"))
col6a5 <- subset(test, subset=seurat_clusters %in% c(3))
test <- merge(col6a5, y = c(col11a1,apoe) , project = "project")
fibcount <- as.matrix(GetAssayData(object = test, slot = "counts"))
write.table(fibcount,"./fibcounts.txt", sep = "\t",quote = F, row.names=T)
saveRDS(col6, file = "./col6.rds")


## 위에서 얻은 fb count, dc count data 합치기
counts <- merge(fibcount, dccount, by=0, all=TRUE) 
counts[is.na(counts)] <- 0  
write.csv(counts,"../expressioncount.csv")
### and then change gene names to ENSEMBLE ID via DAVID GENE ID CONVERSION TOOL and save as counts.txt

## metadata FB, DC : 모든 subcluster 마다 []meta를 만든 후 나중에 전체 merge
Cell <- colnames(fibcount)
col6meta <- as.data.frame(Cell)
cell_type <- c(rep("COL6A5+", length(Cell)))
col6meta$cell_type <- cell_type
met <- rbind(apoemeta,col11meta,col6meta, dc1meta, dc2meta, inflameta)
write.table(met, "../metadata.txt", quote = FALSE, sep = "\t", row.names = F)

## 이후 과정은 linux command
#cellphonedb method statistical_analysis metadata.txt counts.txt --iterations=10 --threads=2
#cellphonedb plot dot_plot
#cellphonedb plot heatmap_plot yourmeta.txt


