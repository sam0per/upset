rm(list = ls())

# PKGs ----
library(ggplot2)
# install.packages("ComplexUpset")
library(ComplexUpset)

# READ ----
mx_fls <- list.files(path = "../../diff_ab_results", pattern = "_int.tsv|_clr.tsv", full.names = TRUE)
mx_fls

mx_dts <- lapply(X = 1:length(mx_fls), FUN = function(x) {
  adt <- read.table(file = mx_fls[x], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
})

# PREP ----
thr <- 0.1

cln <- unique(unlist(lapply(mx_dts, colnames)))
cln <- grep(pattern = "^diagnosis", x = cln, value = TRUE)
(cln <- grep(pattern = "Coefficient$", x = cln, value = TRUE))

mx_sp <- lapply(X = 1:length(mx_fls), FUN = function(x) {
  sdt <- mx_dts[[x]][order(mx_dts[[x]][,cln], decreasing = TRUE),]
  sdt <- sdt[1:(nrow(sdt)*thr),]
  return(sdt)
})

sum(rownames(mx_sp[[1]]) %in% rownames(mx_sp[[2]]))
mx_yn <- data.frame(Taxon = unique(unlist((lapply(mx_sp, rownames)))))

for (i in 1:length(mx_fls)) {
  mx_yn[, basename(mx_fls[i])] <- as.integer(mx_yn$Taxon %in% rownames(mx_sp[[i]]))
}
head(mx_yn)

# UPSET ----
upset(data = mx_yn, intersect = basename(mx_fls), name = "Model", set_sizes = FALSE)
