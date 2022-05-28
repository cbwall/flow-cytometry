## This script applies stochastic neighbor embedding to a FCM classification
## model developed by fcm_model_SG.r or fcm_model_AF.r).  This is a good
## way to test the efficacy of the segmentation model.


# Load in packages
if (!require("pacman")) install.packages("pacman"); library(pacman) # for rapid install if not in library

pacman::p_load("Rtsne", "oce", "ggplot2")


load('SCOOS/SCCOOS_SG.20200704_155059.som.Rdata')

data <- as.data.frame(som.model$data)
deduped.index <- which(duplicated(data) != T)

unit.classif.dedup <- som.model$unit.classif[deduped.index]
data.deduped <- data[deduped.index,]

data.sne <- Rtsne(data.deduped, initial_dims = (dim(data.deduped)[2] - 1),
                  perplexity = 100)

## plot SNE output

flow.col <- oce.colorsFreesurface(k)

point.col <- flow.col[som.cluster[unit.classif.dedup]]

plot(data.sne$Y,
     type = 'n')

points(data.sne$Y,
       col = point.col)

## evaluate correlations between SNE dimensions and parameters

dim.1.cor <- c()
dim.2.cor <- c()

for(var in colnames(data)){
  dim.1.cor <- append(dim.1.cor, cor(data.deduped[,var], data.sne$Y[,1]))
  dim.2.cor <- append(dim.2.cor, cor(data.deduped[,var], data.sne$Y[,2]))
}

dim.cor <- as.data.frame(cbind(dim.1.cor, dim.2.cor))

colnames(dim.cor) <- c('dim.1', 'dim.2')
row.names(dim.cor) <- colnames(data)

