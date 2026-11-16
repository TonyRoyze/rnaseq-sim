# Small base-R helpers shared by the simulation modules.
`%||%` <- function(x, y) if (is.null(x)) y else x

validate_simulation_scalar <- function(x, name, lower = 0, upper = Inf,
                                       include_lower = TRUE) {
  valid_range <- if (include_lower) x >= lower else x > lower
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    valid_range && x <= upper
  if (!ok) {
    range_text <- if (include_lower) paste0("[", lower, ", ", upper, "]") else
      paste0("(", lower, ", ", upper, "]")
    stop(sprintf("%s must be a finite number in %s.", name, range_text))
  }
}
