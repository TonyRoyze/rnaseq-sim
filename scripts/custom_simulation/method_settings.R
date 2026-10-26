# Setting profiles and method adapters for the existing DE method set.
# Each method gets a baseline plus one-factor-at-a-time variants.

make_method_settings <- function(profile = c("core", "quick")) {
  profile <- match.arg(profile)
  settings <- list()
  add <- function(method, id, args = list(), baseline = FALSE, quick = FALSE) {
    settings[[length(settings) + 1L]] <<- list(
      method = method, setting_id = id, args = args,
      baseline = baseline, quick = quick
    )
  }
  add_one_factor <- function(method, base_args, field, values, short_names,
                             quick_first = TRUE) {
    for (i in seq_along(values)) {
      args <- base_args
      args[field] <- list(values[[i]])
      add(method, paste0(field, "-", short_names[i]), args,
          quick = quick_first && i == 1L)
    }
  }

  base <- list(fit_type = "parametric", size_factor = "ratio",
               independent_filtering = TRUE, cooks_cutoff = TRUE)
  add("DESeq2", "default", base, baseline = TRUE)
  add_one_factor("DESeq2", base, "fit_type", c("local", "mean"), c("local", "mean"))
  add_one_factor("DESeq2", base, "size_factor", c("poscounts", "iterate"), c("poscounts", "iterate"))
  add_one_factor("DESeq2", base, "independent_filtering", list(FALSE), "off")
  add_one_factor("DESeq2", base, "cooks_cutoff", list(FALSE), "off")

  base <- list(fit_type = "local")
  add("DESeq", "default", base, baseline = TRUE)
  add_one_factor("DESeq", base, "fit_type", "parametric", "parametric")

  for (method in c("edgeR.lrt", "edgeR.qlf")) {
    base <- list(norm_method = "TMM", robust = FALSE)
    add(method, "default", base, baseline = TRUE)
    add_one_factor(method, base, "norm_method",
                   c("TMMwsp", "RLE", "upperquartile", "none"),
                   c("TMMwsp", "RLE", "upperquartile", "none"))
    add_one_factor(method, base, "robust", list(TRUE), "on")
  }

  base <- list(norm_method = "lr", dispersion_trend = FALSE, equal_var = FALSE)
  add("DSS", "default", base, baseline = TRUE)
  add_one_factor("DSS", base, "norm_method", c("quantile", "total", "median"),
                 c("quantile", "total", "median"))
  add_one_factor("DSS", base, "dispersion_trend", list(TRUE), "trend-on")
  add_one_factor("DSS", base, "equal_var", list(TRUE), "equal-var")

  base <- list(norm_method = "TMM", voom_norm = "none", robust = FALSE, trend = FALSE)
  add("voom", "default", base, baseline = TRUE)
  add_one_factor("voom", base, "norm_method",
                 c("TMMwsp", "RLE", "upperquartile", "none"),
                 c("TMMwsp", "RLE", "upperquartile", "none"))
  add_one_factor("voom", base, "voom_norm", c("quantile"), "quantile")
  add_one_factor("voom", base, "robust", list(TRUE), "robust")
  add_one_factor("voom", base, "trend", list(TRUE), "trend")

  rank_norms <- c("TMM", "TMMwsp", "RLE", "upperquartile")
  base <- list(norm_method = "none", prior_count = 1, exact = NA)
  add("Wilcoxon", "default", base, baseline = TRUE)
  add_one_factor("Wilcoxon", base, "norm_method", rank_norms,
                 c("TMM", "TMMwsp", "RLE", "upperquartile"))
  add_one_factor("Wilcoxon", base, "prior_count", list(0.5), "prior-0.5")
  add_one_factor("Wilcoxon", base, "exact", list(FALSE), "asymptotic")

  base <- list(norm_method = "none", prior_count = 1, equal_var = FALSE)
  add("T.test", "default", base, baseline = TRUE)
  add_one_factor("T.test", base, "norm_method", rank_norms,
                 c("TMM", "TMMwsp", "RLE", "upperquartile"))
  add_one_factor("T.test", base, "prior_count", list(0.5), "prior-0.5")
  add_one_factor("T.test", base, "equal_var", list(TRUE), "equal-var")

  base <- list(norm_method = "AH2010")
  add("NBPSeq", "default", base, baseline = TRUE)
  add_one_factor("NBPSeq", base, "norm_method", list(NULL), "none")
  add_one_factor("NBPSeq", base, "dispersion_model", c("NBP", "NB2"), c("NBP", "NB2"))

  base <- list(norm_method = "tmm", k = 0.5, permutations = 50,
               smoothing = 1.5, filter = 1)
  add("NOISeq", "default", base, baseline = TRUE)
  add_one_factor("NOISeq", base, "norm_method", c("uqua", "n"), c("uqua", "none"))
  add_one_factor("NOISeq", base, "k", list(1), "k-1")
  add_one_factor("NOISeq", base, "permutations", list(100), "perm-100")
  add_one_factor("NOISeq", base, "smoothing", c(0.5, 2), c("smooth-0.5", "smooth-2"))
  add_one_factor("NOISeq", base, "filter", list(0), "filter-off")

  base <- list(norm_method = NULL, replace_outliers = TRUE, use_a_fold = FALSE)
  add("ABSSeq", "package-default", base, baseline = TRUE)
  add_one_factor("ABSSeq", base, "norm_method",
                 c("TMM", "qtotal", "total", "quartile", "geometric"),
                 c("TMM", "qtotal", "total", "quartile", "geometric"))
  add_one_factor("ABSSeq", base, "replace_outliers", list(FALSE), "keep-outliers")
  add_one_factor("ABSSeq", base, "use_a_fold", list(TRUE), "aFold")

  base <- list(norm_method = "TMM", bootstrap = 1000, top_k = NULL)
  add("ROTS", "default", base, baseline = TRUE)
  add_one_factor("ROTS", base, "norm_method", list("none"), "none")
  add_one_factor("ROTS", base, "bootstrap", c(100, 500), c("B-100", "B-500"))
  add_one_factor("ROTS", base, "top_k", list(500), "K-500")

  if (profile == "quick") {
    settings <- Filter(function(x) x$baseline || x$quick, settings)
  }
  settings
}

normalization_factors <- function(counts, method) {
  y <- edgeR::DGEList(counts = counts)
  if (identical(method, "none")) {
    y$samples$norm.factors <- 1
  } else {
    y <- edgeR::calcNormFactors(y, method = method)
  }
  y
}

run_method_setting <- function(dge, setting, seed = 1L) {
  counts <- as.matrix(dge$counts)
  groups <- droplevels(factor(dge$samples$group))
  if (nlevels(groups) != 2L) stop("Each simulated dataset must contain two groups.")
  keep <- rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]
  groups <- groups
  if (!nrow(counts)) stop("No expressed genes remain in this dataset.")
  args <- setting$args
  method <- setting$method
  genes <- rownames(counts)
  p <- NULL

  if (method == "DESeq2") {
    coldata <- data.frame(group = groups, row.names = colnames(counts))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                          design = ~ group)
    dds <- DESeq2::DESeq(dds, fitType = args$fit_type, sfType = args$size_factor,
                         minReplicatesForReplace = if (args$cooks_cutoff) 7 else Inf,
                         quiet = TRUE)
    result <- DESeq2::results(dds, independentFiltering = args$independent_filtering,
                              cooksCutoff = args$cooks_cutoff)
    p <- result$pvalue
  } else if (method == "DESeq") {
    suppressPackageStartupMessages(library(locfit))
    cds <- DESeq::newCountDataSet(counts, groups)
    cds <- DESeq::estimateSizeFactors(cds)
    cds <- DESeq::estimateDispersions(cds, fitType = args$fit_type)
    result <- DESeq::nbinomTest(cds, levels(groups)[1], levels(groups)[2])
    p <- result$pval
    names(p) <- result$id
  } else if (method %in% c("edgeR.lrt", "edgeR.qlf")) {
    y <- normalization_factors(counts, args$norm_method)
    design <- stats::model.matrix(~ 0 + groups)
    fit_disp <- edgeR::estimateDisp(y, design, robust = args$robust)
    if (method == "edgeR.lrt") {
      fit <- edgeR::glmFit(fit_disp, design)
      test <- edgeR::glmLRT(fit, contrast = c(-1, 1))
    } else {
      fit <- edgeR::glmQLFit(fit_disp, design, robust = args$robust)
      test <- edgeR::glmQLFTest(fit, contrast = c(-1, 1))
    }
    p <- edgeR::topTags(test, n = nrow(counts))$table$PValue
    names(p) <- rownames(edgeR::topTags(test, n = nrow(counts))$table)
  } else if (method == "DSS") {
    design_data <- data.frame(
      designs = as.numeric(groups) - 1,
      row.names = colnames(counts)
    )
    seq_data <- DSS::newSeqCountSet(counts, design_data)
    seq_data <- DSS::estNormFactors(seq_data, method = args$norm_method)
    seq_data <- DSS::estDispersion(seq_data, trend = args$dispersion_trend)
    result <- DSS::waldTest(seq_data, 1, 0,
                            equal.var = args$equal_var)
    p <- result$pval
    names(p) <- rownames(result)
  } else if (method == "voom") {
    y <- normalization_factors(counts, args$norm_method)
    design <- stats::model.matrix(~ groups)
    v <- limma::voom(y, design, plot = FALSE,
                     normalize.method = args$voom_norm)
    fit <- limma::eBayes(limma::lmFit(v, design), robust = args$robust,
                         trend = args$trend)
    result <- limma::topTable(fit, coef = ncol(design), n = nrow(counts),
                              sort.by = "none")
    p <- result$P.Value
    names(p) <- rownames(result)
  } else if (method %in% c("Wilcoxon", "T.test")) {
    y <- normalization_factors(counts, args$norm_method)
    effective_lib <- y$samples$lib.size * y$samples$norm.factors
    log_cpm <- log2(sweep(counts, 2, effective_lib, "/") * 1e6 + args$prior_count)
    p <- vapply(seq_len(nrow(log_cpm)), function(i) {
      x <- log_cpm[i, groups == levels(groups)[1]]
      z <- log_cpm[i, groups == levels(groups)[2]]
      if (method == "Wilcoxon") {
        exact <- if (is.na(args$exact)) NULL else args$exact
        stats::wilcox.test(x, z, exact = exact)$p.value
      } else {
        if (length(unique(x)) == 1L && length(unique(z)) == 1L && x[1] == z[1]) return(1)
        stats::t.test(x, z, var.equal = args$equal_var)$p.value
      }
    }, numeric(1))
    names(p) <- genes
  } else if (method == "NBPSeq") {
    norm_factors <- NBPSeq::estimate.norm.factors(
      counts, lib.sizes = colSums(counts), method = args$norm_method
    )
    dispersion_model <- args$dispersion_model
    if (is.null(dispersion_model)) dispersion_model <- "NBQ"
    obj <- NBPSeq::prepare.nbp(counts, groups, lib.sizes = colSums(counts),
                               norm.factors = norm_factors)
    obj <- NBPSeq::estimate.disp(obj, model = dispersion_model, print.level = 0)
    result <- NBPSeq::exact.nb.test(obj, levels(groups)[1], levels(groups)[2],
                                    print.level = 0)
    p <- result$p.values
    names(p) <- genes
  } else if (method == "NOISeq") {
    factors <- data.frame(condition = groups, row.names = colnames(counts))
    input <- NOISeq::readData(data = counts, factors = factors)
    result <- NOISeq::noiseqbio(
      input, norm = args$norm_method, factor = "condition", k = args$k,
      r = args$permutations, adj = args$smoothing, filter = args$filter,
      random.seed = seed, plot = FALSE
    )
    table <- result@results[[1]]
    p <- 1 - table$prob
    names(p) <- rownames(table)
  } else if (method == "ABSSeq") {
    object <- if (is.null(args$norm_method)) {
      ABSSeq::ABSDataSet(counts = counts, groups = groups)
    } else {
      ABSSeq::ABSDataSet(counts = counts, groups = groups,
                         normMethod = args$norm_method)
    }
    object <- ABSSeq::ABSSeq(object, replaceOutliers = args$replace_outliers,
                             useaFold = args$use_a_fold, quiet = TRUE)
    result <- ABSSeq::results(object, c("pvalue"))
    p <- as.numeric(result)
    names(p) <- names(result)
  } else if (method == "ROTS") {
    y <- normalization_factors(counts, args$norm_method)
    expression <- edgeR::cpm(y, log = TRUE, prior.count = 0.5)
    fit <- ROTS::ROTS(expression, groups = as.numeric(groups), B = args$bootstrap,
                      K = args$top_k, seed = seed, log = FALSE, verbose = FALSE)
    result <- fit$data
    p <- fit$pvalue
    names(p) <- rownames(result)
  } else {
    stop("No runner is defined for method: ", method)
  }

  if (is.null(names(p)) || !all(genes %in% names(p))) {
    if (length(p) != length(genes)) stop("Method returned p-values that cannot be aligned to genes.")
    names(p) <- genes
  }
  p <- as.numeric(p[genes])
  p[!is.finite(p) | is.na(p)] <- 1
  names(p) <- genes
  p
}

summarize_method_result <- function(p, truth, alpha = 0.05) {
  truth <- as.logical(truth)
  adjusted <- stats::p.adjust(p, method = "BH")
  called <- adjusted <= alpha
  n_de <- sum(truth)
  n_non_de <- sum(!truth)
  auc <- if (n_de > 0 && n_non_de > 0) {
    ranks <- rank(-p, ties.method = "average")
    (sum(ranks[truth]) - n_de * (n_de + 1) / 2) / (n_de * n_non_de)
  } else NA_real_
  tp <- sum(called & truth)
  fp <- sum(called & !truth)
  data.frame(
    auc = auc,
    true_positive_rate = if (n_de) tp / n_de else NA_real_,
    observed_fdr = if (sum(called)) fp / sum(called) else 0,
    discoveries = sum(called),
    alpha = alpha
  )
}
