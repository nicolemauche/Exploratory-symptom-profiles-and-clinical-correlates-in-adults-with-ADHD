# Supplementary Table S9: proportion of the sample assigned to ED-related classes under the primary
# grouping, posterior weighting, and alternative exploratory groupings.
# Intended for R 4.6.1; requires readxl. No model is re-estimated; only fitted model objects are read.
#
# Inputs (selected interactively):
#   1. IN_EX_Final_WRI_neu.xlsx, sheet IN_EX_Final (28 WRAADDS items, Class_neu = reported six-class assignment)
#   2. fit_k6.rds and 3. fit_k5.rds from the output folder of the LCA workflow (lca/run_final_lca.R)
# Output: folder "alternative_groupings_S9" next to the Excel file. Aggregate counts and proportions only;
#   no participant-level data are written.
#
# Class labels: LC1-LC6 correspond to Class_neu 1-6. The model classes of fit_k6.rds are mapped to Class_neu
# after checking that the item data stored in the model object ($y) equal the Excel items row by row.
# Five-class grouping: the three classes with n = 219, 71 and 198 (unique sizes in the five-class solution).

ED6_CLASSES     <- c(1, 2, 4)                           # primary grouping: LC1, LC2, LC4
SIX_CLASS_SIZES <- c(68, 184, 79, 204, 93, 121)         # LC1-LC6, most likely class assignment
ED5_CLASS_SIZES <- c(219, 71, 198)
groupings <- list(                                      # alternative six-class groupings (LC numbers)
  "Primary six-class grouping" = c(1, 2, 4),
  "Excluding LC4"              = c(1, 2),
  "LC2 only"                   = 2,
  "Including LC6"              = c(1, 2, 4, 6))

message("Running under R ", getRversion(), " (intended: R 4.6.1)")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Please install readxl first: install.packages('readxl')")

message("Select IN_EX_Final_WRI_neu.xlsx");       xlsx <- file.choose()
message("Select fit_k6.rds (LCA output folder)"); f6 <- readRDS(file.choose())
message("Select fit_k5.rds (LCA output folder)"); f5 <- readRDS(file.choose())
x <- as.data.frame(readxl::read_excel(xlsx, sheet = "IN_EX_Final"))
N <- nrow(x)

# ---- checks ----------------------------------------------------------------------
items <- c(paste0("ua", 1:5), paste0("hy", 1:3), paste0("temp", 1:3), paste0("afflab", 1:4),
           paste0("emo", 1:3), paste0("des", 1:5), paste0("imp", 1:5))
if (N != 749) stop("Expected N = 749, found ", N)
if (!all(c(items, "Class_neu") %in% names(x))) stop("Excel file lacks the 28 items or Class_neu.")
X <- as.matrix(x[, items]); dimnames(X) <- NULL
for (nm in c("f6", "f5")) {
  m <- get(nm)
  if (length(m$P) != c(f6 = 6, f5 = 5)[[nm]]) stop(nm, " does not contain the expected number of classes.")
  Y <- as.matrix(m$y); dimnames(Y) <- NULL
  if (!identical(dim(Y), dim(X)) || !all(Y == X))
    stop(nm, ": item data in the model object do not match the Excel items row by row.")
}

# ---- six-class: map model classes to Class_neu (= LC numbering) ------------------
tab6 <- table(model = factor(f6$predclass, levels = 1:6), Class_neu = factor(x$Class_neu, levels = 1:6))
map6 <- apply(tab6, 1, which.max)                       # model class -> LC number
if (sum(apply(tab6, 1, max)) != N || anyDuplicated(map6))
  stop("fit_k6.rds does not reproduce Class_neu one-to-one.")
lc <- map6[f6$predclass]                                # LC number per participant (most likely class)
if (!identical(as.vector(table(factor(lc, levels = 1:6))), as.integer(SIX_CLASS_SIZES)))
  stop("Six-class sizes do not match the reported solution.")
post_lc <- f6$posterior[, order(map6), drop = FALSE]    # posterior columns reordered to LC1-LC6
if (max(abs(rowSums(post_lc) - 1)) > 1e-6) stop("Posterior probabilities do not sum to 1.")

row_modal <- function(label, classes, n) data.frame(
  grouping = label, classes = paste0("LC", classes, collapse = ", "), assignment = "most likely class",
  n = n, percent = round(100 * n / N, 1))

res <- do.call(rbind, lapply(names(groupings), function(g)
  row_modal(g, groupings[[g]], sum(lc %in% groupings[[g]]))))

# posterior-weighted proportion for the primary grouping: per participant, sum of the posterior
# probabilities of LC1, LC2 and LC4; the sum across participants is the expected count
exp_n <- sum(rowSums(post_lc[, ED6_CLASSES, drop = FALSE]))
res <- rbind(res[1, ],
             data.frame(grouping = "Primary six-class grouping, posterior-weighted",
                        classes = paste0("LC", ED6_CLASSES, collapse = ", "), assignment = "posterior-weighted",
                        n = round(exp_n, 1), percent = round(100 * exp_n / N, 1)),
             res[-1, ])

# ---- five-class grouping ---------------------------------------------------------
n5 <- as.vector(table(factor(f5$predclass, levels = 1:5)))
if (anyDuplicated(n5)) stop("Five-class sizes are not unique; classes cannot be identified by size.")
ed5 <- match(ED5_CLASS_SIZES, n5)
if (anyNA(ed5)) stop("Five-class sizes ", paste(n5, collapse = "/"), " do not contain ",
                     paste(ED5_CLASS_SIZES, collapse = "/"), ".")
res <- rbind(res, data.frame(grouping = "Five-class solution",
                             classes = paste0("classes with n = ", paste(ED5_CLASS_SIZES, collapse = ", ")),
                             assignment = "most likely class", n = sum(f5$predclass %in% ed5),
                             percent = round(100 * sum(f5$predclass %in% ed5) / N, 1)))

# emotion-domain means (P of category 3) of the five-class solution, for documentation of the grouping
p3  <- sapply(f5$probs, function(m) m[, 3])
dm5 <- sapply(list(temper = 9:11, affective_lability = 12:15, emotional_overreactivity = 16:18),
              function(j) rowMeans(p3[, j, drop = FALSE]))
dm5 <- data.frame(class5 = 1:5, n = n5, round(dm5, 3), ED_related = 1:5 %in% ed5)


# ---- output --------------------------------------------------------------------------
out <- file.path(dirname(xlsx), "alternative_groupings_S9")
dir.create(out, showWarnings = FALSE)
write.csv(res, file.path(out, "alternative_groupings_S9.csv"), row.names = FALSE)
write.csv(dm5, file.path(out, "five_class_emotion_domain_means.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(res, row.names = FALSE)
cat("\nFive-class emotion-domain means:\n"); print(dm5, row.names = FALSE)
cat("\nOutput written to:", normalizePath(out), "\n")
