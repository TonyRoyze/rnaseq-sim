#!/usr/bin/env Rscript
# Bounded-memory verification for full-factorial v3 simulation output.
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
project_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork=TRUE)), "../.."), mustWork=TRUE)
setwd(project_root)
simulation_version <- "4.0"
opts <- list(output=paste0("results/custom_simulation_v", simulation_version),
  reference="data/DGElist_BRCA.rds", reference_mode="brca",
  references=c("data/DGElist_BRCA.rds","data/DGElist_E-ENAD-34.rds","data/DGElist_GSE150910.rds","data/DGElist_LUSC.rds"),
  replicates=1L, shard_count=1L, batch_size=1000L, n_var=10000L,
  base_seed=20261024L, underdispersion_model="generalized_poisson", de_mode="fixed", de_proportions=c(.1,.2,.3,.4,.5,.6,.75),
  dispersion_estimator="method_of_moments", include_imbalance=TRUE, deep=FALSE)
a <- commandArgs(TRUE)
if (any(a %in% c("--help","-h"))) {
  cat("Usage: Rscript scripts/custom_simulation/check_batch.R [same design flags as run_batch.R] [--deep]\n",
      "Options: --output PATH --reference-mode brca|all|random --reference PATH\n",
      "         --references PATHS --replicates N --shard-count N --n-var N\n",
      "         --underdispersion-model generalized_poisson|com_poisson\n",
      "         --base-seed N --de-mode 1|2|3 --de-proportions LIST\n",
      "         --dispersion-estimator METHOD --include-imbalance|--no-imbalance --batch-size N --deep\n",
      "Default check verifies expected files; --deep also validates each RDS.\n", sep="")
  quit(status=0L)
}
i <- 1L
while (i <= length(a)) {
  k <- a[i]
  if (k == "--deep") { opts$deep <- TRUE; i <- i+1L; next }
  if (k == "--no-imbalance") { opts$include_imbalance <- FALSE; i <- i+1L; next }
  if (k == "--include-imbalance") { opts$include_imbalance <- TRUE; i <- i+1L; next }
  if (!k %in% c("--output","--reference","--reference-mode","--references","--replicates","--shard-count","--batch-size","--n-var","--base-seed","--de-mode","--de-proportions","--dispersion-estimator","--underdispersion-model")) stop("Unknown option: ",k)
  if (i == length(a)) stop("Missing value for ",k)
  v <- a[i+1L]; n <- switch(k,"--output"="output","--reference"="reference","--reference-mode"="reference_mode","--references"="references","--replicates"="replicates","--shard-count"="shard_count","--batch-size"="batch_size","--n-var"="n_var","--base-seed"="base_seed","--de-mode"="de_mode","--de-proportions"="de_proportions","--dispersion-estimator"="dispersion_estimator","--underdispersion-model"="underdispersion_model")
  opts[[n]] <- if(n=="references") trimws(strsplit(v,",",fixed=TRUE)[[1]]) else if(n=="de_proportions") as.numeric(strsplit(v,",",fixed=TRUE)[[1]]) else if(n %in% c("replicates","shard_count","batch_size","n_var","base_seed")) suppressWarnings(as.integer(v)) else v
  i <- i+2L
}
nums <- unlist(opts[c("replicates","shard_count","batch_size","n_var","base_seed")])
if(anyNA(nums)||any(nums<1)) stop("Numeric options must be positive integers.")
if(!opts$underdispersion_model%in%c("generalized_poisson","com_poisson")) stop("Invalid underdispersion model.")
if(!opts$reference_mode%in%c("brca","all","random")) stop("Invalid reference mode.")
if(!opts$dispersion_estimator%in%c("method_of_moments","edgeR")) stop("Invalid dispersion estimator.")
aliases <- c("1"="combination","2"="setting","3"="fixed")
if(opts$de_mode%in%names(aliases)) opts$de_mode <- unname(aliases[opts$de_mode])
if(!opts$de_mode%in%c("combination","setting","fixed")) stop("Invalid DE mode.")
if(!length(opts$de_proportions)||any(!is.finite(opts$de_proportions))||any(opts$de_proportions<=0|opts$de_proportions>=1)||anyDuplicated(opts$de_proportions)) stop("Invalid DE proportions.")
resolve_path <- function(p) if(grepl("^(~|/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p)) path.expand(p) else file.path(project_root,p)
root <- resolve_path(opts$output)
if(!dir.exists(root)) stop("Output directory not found: ",root)
ref_paths <- vapply(if(opts$reference_mode=="brca") opts$reference else opts$references,resolve_path,character(1))
if(any(!file.exists(ref_paths))) stop("Reference files must be present to reproduce job IDs.")
refs <- data.frame(reference_id=tools::file_path_sans_ext(basename(ref_paths)),reference=ref_paths,reference_hash=unname(tools::md5sum(ref_paths)))
source("scripts/custom_simulation/factorial_design.R")
source("scripts/custom_simulation/batch_jobs.R")
source("scripts/custom_simulation/de_proportion_modes.R")
scenarios <- build_scenarios(opts$n_var,opts$include_imbalance,opts$underdispersion_model)
dirs <- list.dirs(root,recursive=FALSE,full.names=TRUE)
dirs <- dirs[grepl("^shard-[0-9]+-of-[0-9]+$",basename(dirs))]
if(!length(dirs)) stop("No shard directories found.")
stamp <- format(Sys.time(),"%Y%m%d-%H%M%S")
missing_file <- file.path(root,paste0("verification-missing-",stamp,".txt"))
invalid_file <- file.path(root,paste0("verification-invalid-",stamp,".txt"))
checked<-missing<-invalid<-interrupted<-0L
n_de <- if(opts$de_mode=="fixed") length(opts$de_proportions) else 1L
n_ref <- if(opts$reference_mode=="all") nrow(refs) else 1L
expected <- as.double(nrow(scenarios))*opts$replicates*n_de*n_ref
cat(sprintf("Checking %.0f expected datasets; deep validation: %s\n",expected,opts$deep)); flush.console()

validate_file <- function(p,j) tryCatch({
  x<-readRDS(p); d<-x$simulated$dge; t<-x$simulated$truth
  parameters <- x$simulated$parameters
  if(!identical(as.character(x$simulation_version),simulation_version)) stop("wrong version")
  if(is.null(d$counts)||!identical(dim(d$counts),c(opts$n_var,sum(c(j$n_group1,j$n_group2))))) stop("wrong matrix dimensions")
  if(!identical(as.integer(table(factor(d$samples$group,levels=levels(d$samples$group)))),as.integer(c(j$n_group1,j$n_group2)))) stop("wrong group sizes")
  if(is.null(t)||nrow(t)!=opts$n_var||sum(t$is_de)!=j$n_diffexp) stop("wrong DE truth")
  if(!identical(as.character(x$scenario$scenario_id[1]),as.character(j$scenario_id))) stop("wrong scenario ID")
  if(!identical(as.integer(x$seed),as.integer(j$seed))) stop("wrong seed")
  if(is.null(parameters)||!identical(as.character(parameters$count_model),as.character(j$count_model))) stop("wrong count model metadata")
  if(!identical(as.integer(parameters$n_diffexp),as.integer(j$n_diffexp))) stop("wrong DE metadata")
  if(!identical(as.integer(parameters$group_sizes),as.integer(c(j$n_group1,j$n_group2)))) stop("wrong group-size metadata")
  if(!identical(as.numeric(parameters$dispersion_scale),as.numeric(j$dispersion_scale))) stop("wrong dispersion metadata")
  if(sum(t$direction=='up')!=round(j$n_diffexp*j$upregulated_fraction)) stop("wrong up/down DE direction counts")
  if(any(t$fold_change_group2_vs_group1[t$direction=='up']<j$fold_change_min) ||
     any(t$fold_change_group2_vs_group1[t$direction=='down']>1/j$fold_change_min) ||
     any(t$fold_change_group2_vs_group1[!t$is_de]!=1)) stop("fold changes do not match DE truth")
  if(any(!is.finite(d$counts))||any(d$counts<0)||any(d$counts!=floor(d$counts))) stop("invalid counts")
  TRUE
},error=function(e) conditionMessage(e))

for(shard in seq_len(opts$shard_count)) {
  opts$shard_index <- shard
  shard_dir <- file.path(root,sprintf("shard-%03d-of-%03d",shard,opts$shard_count))
  it<-make_job_iterator(scenarios,opts,refs)
  cat(sprintf("Shard %d/%d: %.0f expected jobs\n",shard,opts$shard_count,it$assigned)); flush.console()
  repeat {
    batch<-it$next_batch(); if(is.null(batch)) break
    for(q in seq_len(nrow(batch))) {
      j<-batch[q,,drop=FALSE]; nm<-job_output_name(j)
      expected_path <- file.path(shard_dir,nm)
      candidates <- if(file.exists(expected_path)) expected_path else {
        fallback <- file.path(setdiff(dirs,shard_dir),nm)
        fallback[file.exists(fallback)]
      }
      if(!length(candidates)) {
        missing<-missing+1L; cat(expected_path,"\n",file=missing_file,append=TRUE)
        if(any(file.exists(paste0(file.path(dirs,nm),".tmp")))) interrupted<-interrupted+1L
      } else if(opts$deep) {
        ok<-validate_file(candidates[1],j)
        if(!isTRUE(ok)) { invalid<-invalid+1L; cat(candidates[1],"\t",ok,"\n",file=invalid_file,append=TRUE,sep="") }
      }
      checked<-checked+1L
      if(checked%%10000L==0L) {cat(sprintf("  checked %.0f; missing %.0f; invalid %.0f\n",checked,missing,invalid));flush.console()}
    }
    rm(batch)
  }
}
cat(sprintf("\nRESULT expected=%.0f checked=%.0f missing=%d invalid=%d interrupted_tmp=%d\n",expected,checked,missing,invalid,interrupted))
if(missing) cat("Missing paths:",missing_file,"\n")
if(invalid) cat("Invalid datasets:",invalid_file,"\n")
if(!missing&&!invalid) cat("All expected datasets were found and passed the selected checks.\n")
quit(status=if(missing||invalid) 1L else 0L)
