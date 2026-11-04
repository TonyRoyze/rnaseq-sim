# Shared console/file logging. Child processes inherit the parent's log path.
pipeline_log_state <- new.env(parent=emptyenv())
start_pipeline_log <- function(path=Sys.getenv('RNA_PIPELINE_LOG', unset='')) {
  if (!nzchar(path)) return(invisible(NULL))
  con <- file(path, open='at')
  sink(con, split=TRUE)
  pipeline_log_state$con <- con
  options(error=function() {
    cat(format(Sys.time(), '[%Y-%m-%d %H:%M:%S]'), 'ERROR:', geterrmessage(), '\n')
    flush(con)
    quit(save='no', status=1L, runLast=FALSE)
  })
  invisible(con)
}
pipeline_status <- function(stage, text) {
  cat(format(Sys.time(), '[%Y-%m-%d %H:%M:%S]'), stage, text, '\n')
  if (exists("con", envir=pipeline_log_state, inherits=FALSE)) flush(pipeline_log_state$con)
  flush.console()
}
