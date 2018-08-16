#!/usr/bin/env Rscript


# for arandi
library(mcclust) # install.packages("mcclust")

# for %>% 
library(dplyr)

# to call C++ file
library(Rcpp)
sourceCpp("category_cpp.cpp")


phi_prior <- function(matrix_data){
  unnamed_list_prop <- lapply(1:ncol(matrix_data), 
                              function(i) {table(matrix_data[,i])}
                              ) %>% 
    unname %>% 
    lapply('/', nrow(matrix_data))
}

class_priors <- phi_prior(tmp4)

?lapply
qwe <- lapply(1:20, function(i) {table(tmp5[,i])})
x <- unname(qwe)
is.list(x)
qw <- apply(tmp5, 2, table)
y <- unname(qw)
is.list(y)
names

MyData <- read.csv(file="C:/Users/stephen/Desktop/GalactoseData.csv", header=TRUE, sep=",")
MyData <- read.csv(file="C:/Users/stephen/Desktop/pancan12_endometrioid_ovarian.csv", header=TRUE, sep=",")

# BLASPHEMY
MyData -> endometrioid_ovarian_data_discrete

tmp         <- t(endometrioid_ovarian_data_discrete[ -1])
# tmp2        <- cbind(anno_col_endometrioid_ovarian, tmp)
tmp2 <- MyData

# tmp2$cancer <- as.numeric(droplevels(MyData$cancer))-1

myArandiVec <- vector(length = nrow(tmp2))
for(i in 1:nrow(tmp2)){
  myArandiVec[i] <- mcclust::arandi(tmp2$cancer, 
                                    tmp2[,i])
} 

# tmp3 <- tmp2[, sort(myArandiVec, decreasing = T, index.return = T)$ix[2:21]]
# tmp3
pheatmap::pheatmap(tmp3)
tmp3 <- tmp2[,sort(myArandiVec, decreasing = T, index.return = T)$ix[2:21]]
tmp5 <- as.matrix(tmp2[,sort(myArandiVec, decreasing = T, index.return = T)$ix[1:21]])
pheatmap::pheatmap(tmp3)

tmp4 <- as.matrix(tmp3)
class_priors <- phi_prior(tmp4)
p1 <- sum(tmp2$cancer)/nrow(tmp2)
cluster_labels <- sample(c(0,1), nrow(tmp2), replace = T, prob = c(p1, 1- p1))
cluster_weight_priors <- c(0.3, 0.7)
num_clusters <- 2
num_iter <- 10000
burn <- 1000
thinning <- 50
fix_vec <- rep(FALSE, nrow(tmp2))
sim <- sampling(tmp4,
                class_priors,
                cluster_labels,
                fix_vec,
                cluster_weight_priors,
                num_clusters,
                num_iter,
                burn,
                thinning)

dissim <- 1- sim

# Require names to associate data in annotation columns with original data
colnames(dissim) <- MyData$X
rownames(dissim) <- MyData$X

# Example input for annotation_row in pheatmap
annotation_row <- tmp2 %>% dplyr::select(cancer)
rownames(annotation_row) <- rownames(dissim)

# Colour scheme for heatmap
col_pal <- RColorBrewer::brewer.pal(9, "Blues")

# Annotation colours
newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation_row$cancer))))
mycolors <- newCols(length(unique(annotation_row$cancer)))
names(mycolors) <- c(0,1)
mycolors <- list(Class = mycolors)


heat_map <- pheatmap::pheatmap(dissim,
                       annotation_row = annotation_row,
                       annotation_colors = mycolors,
                       cluster_row = T,
                       cluster_cols = T,
                       color = col_pal
)


row_order <- heat_map$tree_row$order
rownames(tmp3) <- rownames(dissim)
heat_data <- tmp3[row_order,]

pheatmap::pheatmap(heat_data, 
                   cluster_rows = F, 
                   cluster_cols = T,
                   annotation_row = annotation_row,
                   annotation_colors = mycolors
                   )
