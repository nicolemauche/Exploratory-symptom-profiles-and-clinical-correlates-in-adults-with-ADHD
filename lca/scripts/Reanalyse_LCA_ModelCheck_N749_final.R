# Model comparison and fit checks for the adult ADHD latent class analysis.
# This script enumerates 1-8 class models and assesses bivariate residuals.
# Run from the package folder in RStudio or from a terminal using this filename.
# Input: IN_EX_Final_WRI_neu.xlsx (sheet IN_EX_Final); output: a user-selected folder.

if (!requireNamespace("readxl", quietly = TRUE)) stop("Install package 'readxl' first.")
if (!requireNamespace("poLCA", quietly = TRUE)) stop("Install package 'poLCA' first.")
if (as.character(getRversion()) != "4.3.3") stop("Requires R 4.3.3; current version is ", getRversion())
if (as.character(utils::packageVersion("poLCA")) != "1.6.0.1") stop("Requires poLCA 1.6.0.1; current version is ", utils::packageVersion("poLCA"))

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "IN_EX_Final_WRI_neu.xlsx"
output_dir <- if (length(args) >= 2) args[[2]] else "lca_modelcheck_n749"
if (!file.exists(input_path) && interactive() && length(args) == 0) {
  message("Select IN_EX_Final_WRI_neu.xlsx")
  input_path <- file.choose()
}
if (!file.exists(input_path)) stop("Workbook not found: ", input_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

items <- c(
  "ua1", "ua2", "ua3", "ua4", "ua5",
  "hy1", "hy2", "hy3",
  "temp1", "temp2", "temp3",
  "afflab1", "afflab2", "afflab3", "afflab4",
  "emo1", "emo2", "emo3",
  "des1", "des2", "des3", "des4", "des5",
  "imp1", "imp2", "imp3", "imp4", "imp5"
)
domains <- c(
  setNames(rep("Inattention", 5), items[1:5]),
  setNames(rep("Hyperactivity", 3), items[6:8]),
  setNames(rep("Temper", 3), items[9:11]),
  setNames(rep("Affective lability", 4), items[12:15]),
  setNames(rep("Emotional overreactivity", 3), items[16:18]),
  setNames(rep("Disorganization", 5), items[19:23]),
  setNames(rep("Impulsivity", 5), items[24:28])
)

raw <- readxl::read_excel(input_path, sheet = "IN_EX_Final")
required <- c("id", "Class_neu", items)
missing <- setdiff(required, names(raw))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
complete <- complete.cases(raw[, items])
dat <- as.data.frame(raw[complete, items])
N <- nrow(dat)
if (N != 749L) stop("Expected N=749 complete cases, found N=", N)
vals <- sort(unique(unlist(dat, use.names = FALSE)))
if (!identical(as.numeric(vals), c(1, 2, 3))) stop("Expected item codes 1, 2, 3.")
if (any(vapply(dat, function(x) any(x != as.integer(x)), logical(1)))) stop("Items must be integer-coded.")

f <- stats::as.formula(paste0("cbind(", paste(items, collapse = ","), ") ~ 1"))
set.seed(123)
nrep_by_k <- c(100L, 100L, 100L, 100L, 100L, 100L, 170L, 130L)
fits <- vector("list", 8)
fit_rows <- vector("list", 8)
for (k in 1:8) {
  message("Fitting ", k, "-class model (", nrep_by_k[k], " random starts)...")
  flush.console()
  fits[[k]] <- poLCA::poLCA(
    f, data = dat, nclass = k, nrep = nrep_by_k[k], maxiter = 100000,
    graphs = FALSE, na.rm = FALSE, verbose = FALSE
  )
  fit <- fits[[k]]
  p <- fit$npar
  aicc <- fit$aic + (2 * p * (p + 1)) / (N - p - 1)
  sabic <- -2 * fit$llik + p * log((N + 2) / 24)
  entropy <- if (k == 1) NA_real_ else {
    1 + sum(fit$posterior * log(pmax(fit$posterior, 1e-12))) / (N * log(k))
  }
  counts <- as.integer(table(factor(fit$predclass, levels = seq_len(k))))
  attempts <- fit$attempts
  n_at_max <- if (is.null(attempts)) NA_integer_ else sum(abs(attempts - fit$llik) <= 0.01, na.rm = TRUE)
  fit_rows[[k]] <- data.frame(
    classes = k, N = N, logLik = fit$llik, parameters = p,
    AIC = fit$aic, AICc = aicc, BIC = fit$bic, SABIC = sabic,
    relative_entropy = entropy,
    min_class_n = min(counts), min_class_percent = 100 * min(counts) / N,
    class_counts = paste(counts, collapse = ","),
    random_starts = if (is.null(attempts)) NA_integer_ else length(attempts),
    starts_at_best_LL = n_at_max,
    percent_starts_at_best_LL = if (is.na(n_at_max)) NA_real_ else 100 * n_at_max / length(attempts)
  )
  # Write each completed model immediately so results survive an interrupted run.
  partial_table <- do.call(rbind, fit_rows[seq_len(k)])
  write.csv(partial_table, file.path(output_dir, "fit_comparison_progress.csv"), row.names = FALSE)
  saveRDS(fits[seq_len(k)], file.path(output_dir, paste0("fits_1to", k, "_final_progress.rds")))
}
fit_table <- do.call(rbind, fit_rows)
write.csv(fit_table, file.path(output_dir, "fit_comparison_1to8.csv"), row.names = FALSE)
saveRDS(fits, file.path(output_dir, "fits_1to8_final.rds"))

# Refit six classes with the final-model setting (500 random starts).
message("Refitting six-class model (500 random starts)...")
set.seed(500)
fit6 <- poLCA::poLCA(
  f, data = dat, nclass = 6, nrep = 500, maxiter = 100000,
  graphs = FALSE, na.rm = FALSE, verbose = FALSE
)
saveRDS(fit6, file.path(output_dir, "fit6_nrep500.rds"))

# Map the arbitrary model class labels to the stored Class_neu labels.
all_permutations <- function(x) {
  if (length(x) == 1L) return(list(x))
  unlist(lapply(seq_along(x), function(i) {
    lapply(all_permutations(x[-i]), function(rest) c(x[i], rest))
  }), recursive = FALSE)
}
ct <- table(
  model_class = factor(fit6$predclass, levels = 1:6),
  reference_class = factor(raw$Class_neu[complete], levels = 1:6)
)
perms <- do.call(rbind, all_permutations(1:6))
match_n <- apply(perms, 1, function(m) sum(ct[cbind(1:6, m)]))
map_to_ref <- perms[which.max(match_n), ]
write.csv(data.frame(model_class=1:6, mapped_Class_neu=map_to_ref,
                     matched_cases=as.integer(ct[cbind(1:6,map_to_ref)])),
          file.path(output_dir,"class6_reference_match.csv"),row.names=FALSE)
if (max(match_n) != N) stop("The 500-start six-class model does not match Class_neu for all 749 cases; inspect class6_reference_match.csv before continuing.")

# Class-conditional probability of category 3, then average within domains.
item_profiles <- do.call(rbind, lapply(items, function(v) {
  pr <- as.matrix(fit6$probs[[v]])
  data.frame(
    model_class = 1:6,
    reference_class = map_to_ref[1:6],
    item = v, domain = unname(domains[v]),
    probability_category_3 = pr[, 3]
  )
}))
write.csv(item_profiles, file.path(output_dir, "class6_item_probabilities.csv"), row.names = FALSE)
domain_profiles <- aggregate(
  probability_category_3 ~ reference_class + domain,
  data = item_profiles, FUN = mean
)
domain_profiles <- domain_profiles[order(domain_profiles$reference_class, domain_profiles$domain), ]
write.csv(domain_profiles, file.path(output_dir, "class6_domain_mean_probabilities.csv"), row.names = FALSE)
emotion_profiles <- subset(domain_profiles, domain %in% c("Temper", "Affective lability", "Emotional overreactivity"))
write.csv(emotion_profiles, file.path(output_dir, "class6_emotion_domain_means.csv"), row.names = FALSE)

# Pairwise bivariate residual diagnostic. Expected pair tables are calculated
# from the fitted class proportions and conditional item probabilities.
pair_result <- function(fit, data, item_a, item_b, K) {
  obs <- as.matrix(table(factor(data[[item_a]], levels = 1:3), factor(data[[item_b]], levels = 1:3)))
  pred <- matrix(0, 3, 3)
  for (cl in seq_len(K)) {
    pa <- as.numeric(fit$probs[[item_a]][cl, 1:3])
    pb <- as.numeric(fit$probs[[item_b]][cl, 1:3])
    pred <- pred + fit$P[cl] * tcrossprod(pa, pb)
  }
  expct <- N * pred
  use <- obs > 0 & expct > 0
  g2 <- 2 * sum(obs[use] * log(obs[use] / expct[use]))
  data.frame(G2 = g2, df = 4, p = stats::pchisq(g2, df = 4, lower.tail = FALSE),
             min_expected = min(expct))
}
all_pairs <- utils::combn(items, 2, simplify = FALSE)
for (k in 1:8) {
  fit <- fits[[k]]
  bvr <- do.call(rbind, lapply(all_pairs, function(pair) {
    r <- pair_result(fit, dat, pair[1], pair[2], k)
    data.frame(
      classes = k, item1 = pair[1], domain1 = unname(domains[pair[1]]),
      item2 = pair[2], domain2 = unname(domains[pair[2]]), r,
      stringsAsFactors = FALSE
    )
  }))
  bvr$p_bonferroni <- pmin(1, bvr$p * nrow(bvr))
  write.csv(bvr, file.path(output_dir, paste0("bivariate_G2_residuals_k", k, ".csv")), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat("\nCompleted. N =", N, "\n")
cat("Class 6 final log-likelihood:", fit6$llik, "\n")
cat("Best one-to-one Class_neu agreement:", max(match_n), "of", N, "\n")
cat("Model comparison, six-class domain means, bivariate residuals, and R version details saved to:\n",
    normalizePath(output_dir), "\n")
