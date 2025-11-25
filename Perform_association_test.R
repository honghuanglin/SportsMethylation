#!/usr/bin/env Rscript

## Association of smoothed methylation proportion with CTE_status
## Usage: Rscript run_cte_methy_glm.R <chr>

suppressPackageStartupMessages(library(data.table))

## ----------------------
## Parse and check args
## ----------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1L || is.na(args[1]) || args[1] == "") {
  stop("Chromosome argument missing. Usage: Rscript run_cte_methy_glm.R <chr>")
}

chr <- as.character(args[1])

## ----------------------
## Input file paths
## ----------------------
pheno_file <- "../Phenotype_CTE.txt"
methy_file <- sprintf("../methy_proportion_smoothed.chr%s.txt", chr)

## ----------------------
## Read data
## ----------------------
pheno <- read.table(
  pheno_file,
  sep     = "\t",
  header  = TRUE,
  quote   = "",
  as.is   = TRUE
)

infile <- fread(
  methy_file,
  sep         = "\t",
  header      = TRUE,
  data.table  = FALSE
)

## ----------------------
## Align phenotype and methylation samples
## ----------------------
idx_samples <- match(pheno$subjid, colnames(infile), nomatch = 0L)
datfile     <- infile[, idx_samples, drop = FALSE]

newpheno <- subset(pheno, subjid %in% colnames(datfile))

## ----------------------
## Run per-CpG logistic regression
## ----------------------
n_cpg   <- nrow(infile)
result  <- matrix(NA_real_, nrow = n_cpg, ncol = 3L)
colnames(result) <- c("beta", "se", "p")

for (i in seq_len(n_cpg)) {

  if (i %% 100 == 0L) {
    message(sprintf("Processing CpG %d / %d at %s", i, n_cpg, Sys.time()))
  }

  cur_cpg <- as.numeric(datfile[i, ])
  cur_cpg <- as.numeric(scale(cur_cpg))  # z-score

  df <- cbind(newpheno, cur_cpg = cur_cpg)
  df <- subset(df, !is.na(cur_cpg))

  fit_sum <- summary(glm(CTE_status ~ cur_cpg,
                         data   = df,
                         family = "binomial"))

  ## Extract coefficient for cur_cpg
  coef_row <- fit_sum$coefficients["cur_cpg", ]
  result[i, ] <- c(
    beta = coef_row["Estimate"],
    se   = coef_row["Std. Error"],
    p    = coef_row["Pr(>|z|)"]
  )
}

## ----------------------
## Prepare output
## ----------------------
start_pos <- infile[, 2] - 1L
end_pos   <- infile[, 2]

out <- data.frame(
  "#chrom" = infile[, 1],
  start    = start_pos,
  end      = end_pos,
  beta     = result[, "beta"],
  se       = result[, "se"],
  p        = result[, "p"],
  row.names = NULL
)

out_file <- sprintf("result.chr%s.txt", chr)

fwrite(
  out,
  file   = out_file,
  sep    = "\t",
  quote  = FALSE
)
