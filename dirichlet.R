

library(mcclust) # install.packages("mcclust")
library(Rcpp)
sourceCpp("category_cpp.cpp")

phi_prior <- function(matrix_data){
  unnamed_list_prop <- matrix_data %>%
    apply(2, table) %>% 
    unname %>% 
    lapply(unname) %>% 
    lapply('/', nrow(matrix_data))
}

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

tmp3 <- tmp2[, sort(myArandiVec, decreasing = T, index.return = T)$ix[2:21]]
tmp3
pheatmap::pheatmap(tmp3)
tmp3 <- tmp2[,sort(myArandiVec, decreasing = T, index.return = T)$ix[1:21]]
pheatmap::pheatmap(tmp3)

tmp3 <- as.matrix(tmp3)
class_priors <- phi_prior(tmp3)
p1 <- sum(tmp2$cancer)/nrow(tmp2)
cluster_labels <- sample(c(0,1), nrow(tmp2), replace = T, prob = c(p1, 1- p1))
cluster_weight_priors <- c(0.3, 0.7)
num_clusters <- 2
num_iter <- 100
burn <- 10
thinning <- 5
fix_vec <- rep(FALSE, nrow(tmp2))
sim <- sampling(tmp3,
                class_priors,
                cluster_labels,
                fix_vec,
                cluster_weight_priors,
                num_clusters,
                num_iter,
                burn,
                thinning)
