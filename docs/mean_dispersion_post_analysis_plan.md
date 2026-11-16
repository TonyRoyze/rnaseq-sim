# Plan: compare mean and dispersion backgrounds in post-analysis

## Objective

Determine whether differential-expression method performance follows the mean-expression profile, the dispersion profile, or both. This is a post-analysis of existing benchmark outputs; it does not require expanding the simulator's full factorial design.

## Implementation plan

1. **Audit whether the comparison is identifiable from saved results.**
   - Check that each reference dataset has results under matching simulation settings.
   - Prefer runs made with `reference-mode all`. In `random` mode, a scenario may have results from only one reference, leaving mean-versus-dispersion contrasts incomplete.
   - Match count model, method, DE proportion, sample sizes, and all other scenario factors.
   - Produce an audit table of references, scenarios, methods, replicates, and missing combinations before calculating performance contrasts.

2. **Summarize reference profiles.**
   - For each reference, calculate gene-level mean and dispersion summaries from the saved parameter cache or original reference data.
   - Record the dispersion estimator, gene filtering rule, reference checksum, and source file.
   - Compare distributions as well as compact summaries so that a small number of highly expressed or variable genes does not dominate a single grand average.

3. **Construct matched background groups.**
   - For references with meaningfully different profiles, classify them as higher/lower mean and higher/lower dispersion using prespecified summary metrics.
   - Where possible, form four groups: low mean/low dispersion, high mean/low dispersion, low mean/high dispersion, and high mean/high dispersion.
   - If the available references do not cover all four groups, use the closest matched pairs and report which contrasts are unavailable or weakly identified.

4. **Use factorial settings as controls.**
   - Compare matched scenarios across reference groups.
   - Begin with a clean subset: no outliers, no extra zeros, no batch effect, balanced groups, and dispersion scale 1.
   - Repeat selected contrasts under stress settings to see whether background effects change with those conditions.
   - Do not average across unrelated settings before checking for setting-specific differences.

5. **Aggregate method performance at the replicate level.**
   - For every method and matched setting, summarize replicate-level AUC and TPR.
   - Also report observed FDR, noting that NOISeq scores are not calibrated p-values and its FDR summaries are descriptive.
   - Preserve failed or missing analyses in the completeness report.
   - Calculate uncertainty across replicates; do not treat gene-level results as independent replicates.

6. **Estimate and report the contrasts.**
   - **Mean-profile contrast:** performance difference between higher- and lower-mean backgrounds at comparable dispersion.
   - **Dispersion-profile contrast:** performance difference between higher- and lower-dispersion backgrounds at comparable mean.
   - **Interaction:** whether the dispersion contrast changes with mean profile, if all four background groups are represented.
   - Plot paired differences by method and metric with uncertainty intervals and the number of matched settings and replicates.

7. **Implement as a separate, reproducible post-analysis.**
   - Add a script under `scripts/de_analysis/` with its R runner in the same task folder that reads analysis RDS files and simulation metadata without modifying the inputs.
   - Write tidy CSV summaries, a comparison plot, and a short README section describing input paths, grouping rules, matching, exclusions, and interpretation limits.

8. **Interpret the comparison cautiously.**
   - Reference datasets can differ in biology, sequencing depth, preprocessing, and other properties in addition to mean and dispersion.
   - Treat results as a comparison of observed reference backgrounds, not as a universal causal decomposition.
   - If available references cannot provide close matches, report that the effects cannot be cleanly separated from the existing runs. A controlled hybrid simulation can then be considered as a follow-up.

## First deliverable

Create a plan-only audit that reports which references, scenario settings, methods, and replicates overlap. Use that audit to decide whether existing outputs support the contrasts before implementing aggregation and plots.
