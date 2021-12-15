rm(list = ls())

# PKGs ----
library(ggplot2)
# install.packages("ComplexUpset")
library(ComplexUpset)
library(grid)
library(optparse)

# OPT ----
option_list = list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "Directory with input files.",
              metavar = "character"),
  make_option(c("-r", "--regex"), type = "character", default = NULL,
              help = "White-spaced keywords for selecting the ordering variable
                (e.g., diagnosis Coefficient; diagnosis FDR_P-value).",
              metavar = "character"),
  make_option(c("-f", "--fig"), type = "character", default = NULL,
              help = "Directory with figures.",
              metavar = "character"),
  make_option(c("-t", "--threshold"), type = "character", default = NULL,
              help = "Float number used for selecting the top rows after variable sorting.
                If negative, sorting is from high to low otherwise it is from low to high.",
              metavar = "character"))

opt_parser = OptionParser(usage = paste("Rscript make_upset.R -d ../..",
                                        "-r diagnosis Coefficient",
                                        "-f ../../plots/upset",
                                        "-r -0.1"),
                          option_list=option_list,
                          description = paste0("Plot upset between different models ",
                                               "using an arbitrary threshold for subsetting input data."))
opt = parse_args(opt_parser)

if (is.null(opt$wdir) | is.null(opt$fig) | is.null(opt$threshold)){
  print_help(opt_parser)
  stop("Missing one or more necessary flags.\n", call.=FALSE)
}

# READ ----
mx_fls <- list.files(path = opt$wdir, full.names = TRUE)

cat("\nThe following files will be intersected ...\n")
basename(mx_fls)

mx_dts <- lapply(X = 1:length(mx_fls), FUN = function(x) {
  adt <- read.table(file = mx_fls[x], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(adt) <- adt$name
  adt$name <- NULL
  return(adt)
})

# PREP ----
# Negative threshold implies sorting from high to low
thr <- opt$threshold
dir_sort <- thr < 0
thr <- abs(thr)

cln <- unique(unlist(lapply(mx_dts, colnames)))
rx <- strsplit(x = "diagnosis Coefficient", split = " ")[[1]]
# re1 <- "diagnosis"
# re2 <- "Coefficient"
cln <- grep(pattern = paste0("^", rx[1]), x = cln, value = TRUE)
(cln <- grep(pattern = paste0(rx[2], "$"), x = cln, value = TRUE))

mx_sp <- lapply(X = 1:length(mx_fls), FUN = function(x) {
  
  if (cln %in% colnames(mx_dts[[x]])) {
    mx_dts[[x]] <- mx_dts[[x]][order(abs(mx_dts[[x]][,cln]), decreasing = dir_sort),]
  }
  sdt <- mx_dts[[x]][1:(nrow(mx_dts[[x]])*thr),]
  return(sdt)
})

sum(rownames(mx_sp[[1]]) %in% rownames(mx_sp[[2]]))
sum(rownames(mx_sp[[1]]) %in% rownames(mx_sp[[8]]))
mx_yn <- data.frame(Taxon = unique(unlist((lapply(mx_sp, rownames)))))

for (i in 1:length(mx_fls)) {
  mx_yn[, basename(mx_fls[i])] <- as.integer(mx_yn$Taxon %in% rownames(mx_sp[[i]]))
}
head(mx_yn)

# UPSET ----
mx_up <- upset(data = mx_yn, intersect = basename(mx_fls), name = "Model", set_sizes = FALSE) +
  labs(caption = paste0("Mixed models: Top ", thr*100, "% with largest ", rx[1], " ", rx[2],
                        "\nML: Top ranked ", thr*100, "%"))
mx_up
# grid.text("My_Title",x = 0.65, y=0.95, gp=gpar(fontsize=20))

up_fig <- opt$fig
dir.create(up_fig)
ggsave(filename = paste0(up_fig, "/upset_ML_mixed_models_", re1, "_", re2, "_top", thr*100, "perc.pdf"),
       plot = mx_up, dpi = "print")
