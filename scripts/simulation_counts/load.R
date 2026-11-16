# Public entry point for the configurable simulator.
# Source this file from the project root; it links the readable modules in
# scripts/simulation_counts/ and makes their functions available.
module_dir <- dirname(normalizePath(sys.frame(1)$ofile))
for (module in c(
  "helpers.R",
  "estimate_parameters.R",
  "prepare_design.R",
  "generalized_poisson.R",
  "generate_counts.R",
  "simulate_counts.R"
)) {
  source(file.path(module_dir, module), local = FALSE)
}
rm(module_dir, module)
