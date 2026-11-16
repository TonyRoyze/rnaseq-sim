source('scripts/simulation_pipeline/machine_config.R')
for(profile in c('quick','core')) {
  all <- make_method_settings(profile)
  for(method in unique(vapply(all,`[[`,character(1),'method'))) {
    selected <- select_method_settings(all,method)
    expected <- all[vapply(all,function(spec) spec$method==method,logical(1))]
    stopifnot(identical(selected,expected))
  }
  selected <- select_method_settings(all,c('NBPSeq','DESeq2'))
  stopifnot(setequal(unique(vapply(selected,`[[`,character(1),'method')),c('NBPSeq','DESeq2')))
}
for(methods in list(character(),'',c('NBPSeq','NBPSeq'),'nbpseq','UNKNOWN'))
  stopifnot(inherits(try(select_method_settings(all,methods),silent=TRUE),'try-error'))
base <- readLines('config/machine.dcf')
base <- base[seq_len(which(base == 'Methods:')[1]-1L)]
fixture <- c(base, 'Methods:', '    NBPSeq', '    # DESeq2')
lines <- fixture
path <- tempfile(fileext='.dcf')
writeLines(lines,path)
stopifnot(identical(read_machine_config(path)$Methods,'NBPSeq'))
lines <- sub('    # DESeq2','    DESeq2',lines,fixed=TRUE)
writeLines(lines,path)
stopifnot(identical(read_machine_config(path)$Methods,c('NBPSeq','DESeq2')))
writeLines(sub('    NBPSeq','    # NBPSeq',fixture,fixed=TRUE),path)
stopifnot(inherits(try(read_machine_config(path),silent=TRUE),'try-error'))
unlink(path)
cat('PASS: uncomment selection, multi-method pools, full profile preservation and invalid/empty selections\n')

core <- make_method_settings('core')
expected <- c(ABSSeq=6L, DESeq=1L, DESeq2=3L, DSS=4L, NBPSeq=2L,
              NOISeq=4L, ROTS=2L, T.test=5L, Wilcoxon=5L,
              edgeR.lrt=5L, edgeR.qlf=5L, voom=5L)
observed <- table(vapply(core, `[[`, character(1), 'method'))
stopifnot(identical(as.integer(observed[names(expected)]), unname(expected)))
for (spec in core) if (!spec$baseline)
  stopifnot(grepl('^(norm_method|size_factor)-', spec$setting_id))
cat('PASS: exactly 47 requested method/normalization labels, no other parameter variants\n')
