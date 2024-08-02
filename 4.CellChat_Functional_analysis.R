
9. INTERACTION AND GENE ONTOLOGY 

# Cellchat chino (la mayoría es su script pero fue necesario adaptar alguna cosa)
BiocManager::install("ComplexHeatmap")
library (ComplexHeatmap) 
devtools::install_github("sqjin/CellChat")
library(CellChat)
library(Seurat)
library(reticulate)

#Part I: Create a CellChat object
# El objeto debe de tener una layer de data pero tener anotación
# 1. Create a CellChat object from a Seurat RDS file
View(SCT_objCCA@meta.data)
data.input <- GetAssayData(SCT_objCCA, assay = "SCT", layer = "data") 
Idents(SCT_objCCA)  <- "singler_labels"
labels <- Idents(SCT_objCCA)
meta <- data.frame(group = labels, row.names = names(labels)) 
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
# [1] "Create a CellChat object from a data matrix"
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  Neurons Oligodendrocytes Astrocytes Neurons activated Microglia aNSCs Fibroblasts Fibroblasts activated Monocytes Cardiomyocytes Dendritic cells qNSCs Astrocytes activated Endothelial cells Ependymal Adipocytes NPCs B cells Macrophages Fibroblasts senescent Macrophages activated Microglia activated NK cells T cells 



# Part II: Set and Update CellChatDatabases 
# Load the ligand-receptor interaction database
CellChatDB.human <- CellChatDB.human # Load the human database 
CellChatDB.mouse <- CellChatDB.mouse # Load the mouse database

# Show the ligand-receptor categories
showDatabaseCategory(CellChatDB.human)
showDatabaseCategory(CellChatDB.mouse)

# Set all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB.mouse

# Set a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling")
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "ECM-Receptor")
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Cell-Cell Contact")

# Update CellChatDB
interaction <- CellChatDB.mouse$interaction
complex <- CellChatDB.mouse$complex
cofactor <- CellChatDB.mouse$cofactor
geneInfo <- CellChatDB.mouse$geneInfo

write.csv(interaction,  file = "C:/Users/plaza/OneDrive/Escritorio/TFM/R course/interaction.csv")
write.csv(complex,  file = "C:/Users/plaza/OneDrive/Escritorio/TFM/R course/complex.csv")
write.csv(cofactor,  file = "C:/Users/plaza/OneDrive/Escritorio/TFM/R course/cofactor.csv")
write.csv(geneInfo,  file = "C:/Users/plaza/OneDrive/Escritorio/TFM/R course/geneInfo.csv")

# Adding users’ curated ligand-receptor pairs
# We can add to the csv we have created new things to use for our data

# Update
options(stringsAsFactors = FALSE)
interaction <- read.csv(file = 'C:/Users/plaza/OneDrive/Escritorio/TFM/R course/interaction.csv', 
                        row.names = 1)
complex <- read.csv(file = 'C:/Users/plaza/OneDrive/Escritorio/TFM/R course/complex.csv', 
                    row.names = 1)
cofactor <- read.csv(file = 'C:/Users/plaza/OneDrive/Escritorio/TFM/R course/cofactor.csv', 
                     row.names = 1)
geneInfo <- read.csv(file = 'C:/Users/plaza/OneDrive/Escritorio/TFM/R course/geneInfo.csv', 
                     row.names = 1)

CellChatDB.mouse.updated <- list()
CellChatDB.mouse.updated$interaction <- interaction
CellChatDB.mouse.updated$complex <- complex
CellChatDB.mouse.updated$cofactor <- cofactor
CellChatDB.mouse.updated$geneInfo <- geneInfo

CellChatDB.use <- CellChatDB.mouse.updated


# Part III: Inference of cell-cell communication network 
library(CellChat)
library(patchwork)
library(circlize)
options(stringsAsFactors = FALSE)

# 1.Add CellChatDB in your cellchat object
# Set the full human ligand-receptor interaction database
CellChatDB.use <- CellChatDB.mouse.updated

# Using a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB.mouse.updated, search = "Secreted Signaling")

# Add the Secreted Signaling database in the CellChat object
cellchat@DB <- CellChatDB.use

# 2.Subset and pre-processing the expression data 
# subset the expression data to use less RAM
cellchat <- subsetData(cellchat)

# Pre-processing the expression data
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Optional: project gene expression data onto protein-protein interaction (PPI)
cellchat <- projectData(cellchat, PPI.human)

# 3. Compute the communication probability and infer cellular communication network
# cellchat <- computeCommunProb(cellchat) (if we did not did the prior step)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE) # use the projected data

# 4. Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  aNSCs Fibroblasts Monocytes Cardiomyocytes Dendritic cells Adipocytes B cells Macrophages Fibroblasts senescent Macrophages activated NK cells T cells 

# 5. Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# 6. Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight

# 7. visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()

#  examine the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(2, 6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

mat <- cellchat@net$weight
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}




# Part IV: Visualization of cell-cell communication network
cellchat@netP[["pathways"]]

extractEnrichedLR(cellchat, signaling = c(cellchat@netP[["pathways"]]),
                  geneLR.return = TRUE)
# $pairLR
# interaction_name
# 1     PSAP_GPR37L1
# 2       NRG1_ERBB4
# 3       NRG3_ERBB4
# 4       PTN_PTPRZ1
# 5         PTN_SDC2
# 6         PTN_SDC4
# 7       IGF1_IGF1R
# 8       IL34_CSF1R
# 
# $geneLR
# [1] "Psap"    "Gpr37l1" "Nrg1"    "Nrg3"    "Erbb4"   "Ptn"     "Ptprz1"  "Sdc2"    "Sdc4"    "Igf1"   
# [11] "Igf1r"   "Il34"    "Csf1r"  

# visualize the contribution of each LR pairs to the communication network
netAnalysis_contribution(cellchat, 
                         signaling = c(cellchat@netP[["pathways"]]), 
                         title = "Contribution of each LR pairs")

netAnalysis_contribution(cellchat, 
                         signaling = c(cellchat@netP[["pathways"]][1:5]))

# See the contribution of individual pathways
extractEnrichedLR(cellchat, signaling = "TNF", geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = "TNF")

extractEnrichedLR(cellchat, signaling = "CCL", geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = "CCL")

# Circle plot
netVisual_aggregate(cellchat, signaling = "CCL", layout = "circle")
netVisual_individual(cellchat, signaling = "CCL", layout = "circle")
netVisual_individual(cellchat, signaling = "CCL", 
                     pairLR.use = "CCL5_CCR1",
                     layout = "circle")

# Chord diagram
par(mfrow = c(1, 1), xpd=TRUE)
par(cex = 0.5)
netVisual_aggregate(cellchat, signaling = "CCL", layout = "chord")
netVisual_chord_cell (cellchat, signaling = "CCL")
netVisual_chord_gene (cellchat, signaling = "CCL")

# Chord diagram: group cell clusters into fibroblast, DC and TC cells 
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
names(group.cellType) <- levels(cellchat@idents)
par(mfrow = c(1, 1), xpd=TRUE)
par(cex = 0.5)
netVisual_chord_cell(cellchat, signaling = "CCL", 
                     group = group.cellType, 
                     title.name = paste0("CCL_", "signaling network"))

# Chord diagram: define source and target cell types
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:6), 
                     lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8,
                     lab.cex = 0.5, legend.pos.x = 15)

# Chord diagram: show LR pairs associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8,
                     signaling = c("CCL","CXCL"),legend.pos.x = 8)

# Hierarchy plot 
vertex.receiver = seq(1,4) # define the left portion of cell groups  
netVisual_aggregate(cellchat, signaling = "CCL", 
                    vertex.receiver = vertex.receiver, layout = "hierarchy")
netVisual_individual(cellchat, signaling = "CXCL", 
                     pairLR.use = "CXCL12_CXCR4", 
                     vertex.receiver = vertex.receiver, 
                     layout = "hierarchy")

# Heatmap
netVisual_heatmap(cellchat, signaling = "CCL", color.heatmap = "Reds")

# Violin plot 
plotGeneExpression(cellchat, signaling = "CCL")

# bubble plot 
# bubble plot: show all LR pairs from source to target cell groups
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:12), 
                 remove.isolate = FALSE) 

# bubble plot: show LR pairs associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:12), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)

# Part V: Systematic analysis of cell-cell communication networks
library(CellChat)
library(NMF)
library(ggalluvial)

# 1. Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Scatter plot to visualize aggregated communication networks for each cell type
netAnalysis_signalingRole_scatter(cellchat) # all signaling pathways

# Scatter plot to Visualize selected communication networks
netAnalysis_signalingRole_scatter(cellchat, signaling = "TNF")
netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))

# Heatmap to visualize dominant cell types for each signaling pathway
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 11)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 11)

# Visualize selected outgoing/incoming signals and contributing cell types
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                  signaling = c("CXCL", "CCL"))
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                  signaling = c("CXCL", "CCL"))

# Heatmap to visualize major signaling roles of different cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "CCL", width = 10, 
                                  height = 5, font.size = 10)

# 2. Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate

# Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat, pattern = "outgoing") # infer the number of patterns, NMF
nPatterns = 3 # a suitable number of patterns is the one begin to drop suddenly.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns, width = 5, height = 9)

netAnalysis_river(cellchat, pattern = "outgoing") # river plot
netAnalysis_dot(cellchat, pattern = "outgoing") # dot plot

## Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat,pattern = "incoming", 
                                          k = nPatterns, width = 5, height = 9)

netAnalysis_river(cellchat, pattern = "incoming") # river plot
netAnalysis_dot(cellchat, pattern = "incoming") # dot plot

# EL ANALISIS FUNCIONAL Y ESTRUCTURAL REQUIERE PYTHON. SALTA ERROR
# # 3. Groups signaling pathways based on their functional/structural similarities
# # Identify signaling groups based on functional similarity
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering (cellchat, type = "functional", do.parallel = FALSE)
# 
# # # Visualization in 2D-space
# netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
# 
# # Identify signaling groups based on structure similarity multimeric ligand-receptor complexes, soluble agonists and antagonists, stimulatory and inhibitory co-ligands and co-receptors
# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE)
# 
# # Visualization in 2D-space
# netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


# Part VI: Compare Cell-cell Communication Networks across Biological Conditions
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# 1. Create the NL CellChat object 
SCT_objCCA
table(SCT_objCCA[["sample"]])
data.input = SCT_objCCA@assays$SCT$data # normalized data matrix
meta = SCT_objCCA@meta.data # a dataframe with rownames containing cell meta data
cell.use = rownames(meta)[meta$sample == "ctrl"] #we are only looking into the 'normal' cells
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
unique(meta$singler_labels) # check the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singler_labels")

# 2. Add the Secreted Signaling database in the CellChat object
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# 3.Subset and pre-processing the expression data 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 4. project gene expression data onto protein-protein interaction (PPI)
cellchat <- projectData(cellchat, PPI.human) # PPI.mouse for mouse samples

# 5. Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = FALSE) 

# 6. Filter out the cell-cell communication min.cells = 10
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 7. Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# 8. Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# 9. Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 10. Identify and visualize outgoing communication pattern of secreting cells
selectK(cellchat, pattern = "outgoing") 
nPatterns = 2 
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns, width = 5, height = 9)

# 11. Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat,pattern = "incoming", 
                                          k = nPatterns, width = 5, height = 9)

# Línea 274 problema
# # 12. Identify signaling groups based on functional similarity
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering (cellchat, type = "functional", do.parallel = FALSE)
# 
# # 13. Identify signaling groups based on structure similarity
# cellchat <- computeNetSimilarity(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE)

# 14. Save ctrl cellchat object
saveRDS(cellchat, file="...")

# 15. load both NL and LS objects. La parte de arriba solo te analiza el ctrl, entiendo que deberías de hacer lo mismo con el stim y una vez analizados mergeas los dos cellchat generados
cellchat.ctrl <- readRDS("...")
cellchat.stim <- readRDS("...")
object.list <- list(ctrl = cellchat.ctrl,stim = cellchat.stim)
names(object.list)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat




# Part VII: Continuation of Part VI. Using the merged cellchat Compare Communication Networks across Biological Conditions

m(cellchat.ctrl, cellchat.stim)

# 1. Compare the overall information flow of each signaling pathway
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) # 42 & 47
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

# 2. Compare the total number of interactions and interaction strength
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

###### scatter plot 
# 1.Compare outgoing/incoming interaction strength for all the cell types
count.sum <- sapply(object.list, function(x) {rowSums(x@net$count) + 
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(count.sum), max(count.sum)) # control the dot size 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

# 2. identify signalling changes associated with one cell group 
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC", 
                                     signaling.exclude = "MIF")

###### Circle plots
# 1. show the number of interactions between any two cell populations 
# compute the maximum number of cells and the maximum number of interactions 
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F,
                   edge.weight.max = weight.max[2], edge.width.max = 12, arrow.size = 0.1,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# 2. selected pathway
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), 
                           attribute =c("CXCL"))

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = c("CXCL"), layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10, arrow.size = 0.05, 
                      signaling.name = paste("CXCL", names(object.list)[i]))
}

# 3. Show differential number of interactions or interaction strength among 
# different cell populations, red(increased signaling)/blue(decreased signaling)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), measure = "count", 
                          weight.scale = T, arrow.size = 0.1)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), measure = "weight", 
                          weight.scale = T, arrow.size = 0.1)

# 4. simplify the complicated network to the cell type level
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {
  mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), 
                           attribute = c("idents", "count", "count.merged"))

# show the number of interactions or interaction strength.
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, 
                   label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, 
                   arrow.size = 0.1,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(1, 2),
                          arrow.size = 0.1, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(1, 2),
                          arrow.size = 0.1, measure = "weight.merged", label.edge = T)

# Part VIII: Continuation of Parts VI & VII
###### Heatmap
# 1. Compare outgoing/incoming signaling associated with each cell population
# combining all the identified signaling pathways from different datasets 
all_pathways <- union(object.list[[1]]@netP$pathways, 
                      object.list[[2]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", 
                                        signaling = all_pathways, title = names(object.list)[1],  
                                        width = 5, height = 11, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", 
                                        signaling = all_pathways, title = names(object.list)[2], 
                                        width = 5, height = 11, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht3 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern ="outgoing", 
                                        signaling = all_pathways, title = names(object.list)[1], 
                                        width = 5, height = 11)
ht4 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern ="outgoing", 
                                        signaling = all_pathways, title = names(object.list)[2], 
                                        width = 5, height = 11)
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

ht5 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", 
                                        signaling = all_pathways, title = names(object.list)[1], 
                                        width = 5, height = 11, color.heatmap = "GnBu")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern ="incoming",
                                        signaling = all_pathways, title = names(object.list)[2], 
                                        width = 5, height = 11, color.heatmap = "GnBu")
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))

# 2. selected pathways
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = c("CXCL"), 
                               title.name = paste("CXCL", "signaling ",names(object.list)[i]),
                               color.heatmap = "Reds")
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# 3. show differential interaction number & interaction strength using heatmap
gg1 <- netVisual_heatmap(cellchat, comparison = c(1, 2), measure = "count")
gg2 <- netVisual_heatmap(cellchat, comparison = c(1, 2), measure = "weight")
gg1 + gg2

###### Bubble plots
# 1. compare communication probabilities mediated by ligand-receptor pairs from 
# selected sources and targets cell groups 
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(6:11),  
                 comparison = c(1, 2), angle.x = 45)

# 2. identify the up-regulated ligand-receptor pairs 
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(6:11),  
                 comparison = c(1, 2), max.dataset = 2, 
                 title.name = "Increased signaling in LS", angle.x = 45, 
                 remove.isolate = T)

# 3. identify down-regulated ligand-receptor pairs 
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(6:11),  
                 comparison = c(1, 2), max.dataset = 1, 
                 title.name = "Decreased signaling in LS", angle.x = 45, 
                 remove.isolate = T) 

# Violin plot
View(cellchat@meta)

plotGeneExpression(cellchat, signaling = "CXCL", split.by = "condition", 
                   colors.ggplot = T)




# Part IX: Analyse Cell-cell Communication Networks for Sequencing Based Spatial Dataset (no nos interesa)
