mx_n <- unique(sapply(mx_dts, nrow))

if (length(mx_n)>1) {
  mx_n <- max(mx_n)
}

mx_sp <- matrix(nrow = round(thr*mx_n), ncol = length(mx_fls))
colnames(mx_sp) <- c("Mixed CLR", "Mixed INT")

strsplit(x = gsub(pattern = ".tsv", replacement = "", x = basename(mx_fls)), split = "_")