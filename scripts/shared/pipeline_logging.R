# Shared console/file logging. Child processes inherit the parent's log path.
pipeline_log_state <- new.env(parent=emptyenv())
start_pipeline_log <- function(path=Sys.getenv('RNA_PIPELINE_LOG', unset='')) {
  con <- NULL
  if (nzchar(path)) {
    con <- file(path, open='at')
    sink(con, split=TRUE)
    pipeline_log_state$con <- con
  }
  options(error=function() {
    pipeline_diagnostic('FATAL', geterrmessage(), sys.calls())
    quit(save='no', status=1L, runLast=FALSE)
  })
  invisible(con)
}
pipeline_context <- function(text) {
  pipeline_log_state$context <- text
}
pipeline_diagnostic <- function(level, message, calls) {
  context <- pipeline_log_state$context
  frames <- vapply(tail(calls, 12L), function(call) {
    paste(deparse(call, width.cutoff=100L, nlines=1L), collapse=' ')
  }, character(1))
  lines <- c(sprintf('[%s] %s pid=%d%s | %s',
                     format(Sys.time(), '%Y-%m-%d %H:%M:%S %z'), level, Sys.getpid(),
                     if (is.null(context)) '' else paste0(' | ', context),
                     trimws(message)),
             paste0('  at ', frames))
  writeLines(lines, stderr())
  flush(stderr())
  con <- pipeline_log_state$con
  if (!is.null(con)) {
    writeLines(lines, con)
    flush(con)
  }
}
log_pipeline_errors <- function(expr) {
  withCallingHandlers(expr, error=function(e) {
    pipeline_diagnostic('R ERROR', conditionMessage(e), sys.calls())
  })
}
run_logged_process <- function(executable, args, log, label) {
  pipeline_status('PROCESS START', sprintf('%s | launcher pid=%d | log=%s',
                                           label, Sys.getpid(), log))
  started <- proc.time()[['elapsed']]
  code <- tryCatch(system2(executable, args, stdout=log, stderr=log), error=identity)
  outcome <- if (inherits(code, 'error')) {
    paste('launch error:', conditionMessage(code))
  } else {
    paste('exit status', code)
  }
  detail <- sprintf('%s | launcher pid=%d | %s | elapsed %.3f sec | log=%s',
                    label, Sys.getpid(), outcome, proc.time()[['elapsed']] - started, log)
  cat(sprintf('\n[%s] PROCESS EXIT %s\n', format(Sys.time(), '%Y-%m-%d %H:%M:%S %z'),
              detail), file=log, append=TRUE)
  pipeline_status('PROCESS EXIT', detail)
  if (inherits(code, 'error') || code != 0L) stop(detail, call.=FALSE)
  invisible(code)
}
pipeline_status <- function(stage, text) {
  cat(format(Sys.time(), '[%Y-%m-%d %H:%M:%S]'), stage, text, '\n')
  if (exists("con", envir=pipeline_log_state, inherits=FALSE)) flush(pipeline_log_state$con)
  flush(stdout())
  flush.console()
}
