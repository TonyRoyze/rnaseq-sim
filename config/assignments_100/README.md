# 100 reusable single-method assignments

These are 100 work assignments, not 100 physical machines. Every config runs
one DE method with all its requested core normalization variants. It keeps
30 replicates, 10,000 genes, six local workers, and seeded random selection
among the same four reference files. Each method covers dataset IDs 1–136080
exactly once. Different methods deliberately regenerate the same datasets.

The 27/27/10 assignments for NBPSeq/DESeq2/DESeq plus 36 assignments divided
equally among the other nine methods extend the provisional 8/8/3/11 physical
machine plan. This allocation is not a measured throughput optimum; job counts
are not runtime estimates. Each assignment has a unique Output directory.

| Configs | Method | Assignments | Datasets per assignment |
| --- | --- | --- | --- |
| 001–027 | NBPSeq | 27 | 5,040 |
| 028–054 | DESeq2 | 27 | 5,040 |
| 055–064 | DESeq | 10 | 13,608 |
| 065–068 | edgeR.lrt | 4 | 34,020 |
| 069–072 | edgeR.qlf | 4 | 34,020 |
| 073–076 | DSS | 4 | 34,020 |
| 077–080 | voom | 4 | 34,020 |
| 081–084 | Wilcoxon | 4 | 34,020 |
| 085–088 | T.test | 4 | 34,020 |
| 089–092 | NOISeq | 4 | 34,020 |
| 093–096 | ABSSeq | 4 | 34,020 |
| 097–100 | ROTS | 4 | 34,020 |

## Initial assignment to your 30 machines

| Physical machine(s) | Initial config(s) | Method |
| --- | --- | --- |
| 01–08 | 001–008, respectively | NBPSeq |
| 09–16 | 028–035, respectively | DESeq2 |
| 17–19 | 055–057, respectively | DESeq |
| 20 | 065 | edgeR.lrt |
| 21 | 069 | edgeR.qlf |
| 22 | 073 | DSS |
| 23 | 077 | voom |
| 24 | 081 | Wilcoxon |
| 25 | 085 | T.test |
| 26–27 | 089–090, respectively | NOISeq |
| 28 | 093 | ABSSeq |
| 29–30 | 097–098, respectively | ROTS |

When a machine finishes, give it the next unassigned config for that method.
For example, a free NBPSeq machine takes 009 onward, a free DESeq2 machine
036 onward, and the edgeR.lrt machine takes 066 after completing 065.
The two extra machines in the remaining-method group initially serve NOISeq
and ROTS; change this allocation after measuring representative runtimes.
When a method's assignments finish, a free machine can take another method's
unassigned config if its dependencies are installed.

A single coordinator should record each assignment before starting it.
Separate local output folders on two computers cannot detect duplicate work.
The default `config/machine.dcf` is unchanged; use `--config` explicitly.

## Run a config (Windows PowerShell, from the repository root)

```powershell
.\automations\run_machine.ps1 --config config/assignments_100/config-001.dcf --plan-only
.\automations\run_machine.ps1 --config config/assignments_100/config-001.dcf
$LASTEXITCODE
```

Use the assigned three-digit config number. Physical machine 09 initially
runs config 028, rather than config 009. On macOS/Linux use
`./automations/run_machine.sh` with the same arguments.

Wait for the command to return. Mark Done only after exit code 0 and every
worker reports complete; a silent console alone does not mean completion.
To retry after failure, rerun the same config on the same computer with its
existing Output folder; saved results are reused. Preserve each completed
output folder when starting the next assignment. To move an interrupted
assignment, stop its old workers and transfer its entire Output folder to
the replacement computer, or explicitly choose to recompute it.
Do not change job bounds, methods, seeds, scientific settings or Output after
starting a config. Workers can be reduced before first use on a low-RAM machine;
changing a saved configuration is rejected by the resume snapshot.

Every computer needs the same code/package versions and the four reference
RDS files. Six workers are a starting setting, subject to CPU/RAM capacity.
All datasets and method results are retained; provide adequate disk space and
collect the uniquely named output folders centrally. This package adds no
automatic dispatch, central queue or cross-machine lock.

## Method loading and process overhead

Analysis jobs are filtered to the selected method before execution. Package
namespaces are invoked inside the selected runner branch, rather than loading
all DE packages upfront. Shared dependencies still load (for example edgeR
for count objects and normalization, and limma for voom). A single-method config
does not mean a single R package. The machine launcher passes `--no-restore`;
prepare the pinned package library once before running assignments.

Each worker still starts fresh generation and analysis R processes for each
20-dataset batch. Package loading is reused within an analysis process, but
repeats in subsequent batches. These configs do not implement persistent
analysis processes or additional logging.

## Assignment register

Keep one authoritative coordinator's copy of this register. Status can be
Unassigned, Running, Failed or Done. Do not edit independent copies on all
computers and assume they synchronize. Config numbers are assignment IDs,
not dataset IDs or physical machine IDs.

| Config | Method | Job-Start | Job-End | Assigned machine | Status |
| --- | --- | --- | --- | --- | --- |
| [001](config-001.dcf) | NBPSeq | 1 | 5040 | | Unassigned |
| [002](config-002.dcf) | NBPSeq | 5041 | 10080 | | Unassigned |
| [003](config-003.dcf) | NBPSeq | 10081 | 15120 | | Unassigned |
| [004](config-004.dcf) | NBPSeq | 15121 | 20160 | | Unassigned |
| [005](config-005.dcf) | NBPSeq | 20161 | 25200 | | Unassigned |
| [006](config-006.dcf) | NBPSeq | 25201 | 30240 | | Unassigned |
| [007](config-007.dcf) | NBPSeq | 30241 | 35280 | | Unassigned |
| [008](config-008.dcf) | NBPSeq | 35281 | 40320 | | Unassigned |
| [009](config-009.dcf) | NBPSeq | 40321 | 45360 | | Unassigned |
| [010](config-010.dcf) | NBPSeq | 45361 | 50400 | | Unassigned |
| [011](config-011.dcf) | NBPSeq | 50401 | 55440 | | Unassigned |
| [012](config-012.dcf) | NBPSeq | 55441 | 60480 | | Unassigned |
| [013](config-013.dcf) | NBPSeq | 60481 | 65520 | | Unassigned |
| [014](config-014.dcf) | NBPSeq | 65521 | 70560 | | Unassigned |
| [015](config-015.dcf) | NBPSeq | 70561 | 75600 | | Unassigned |
| [016](config-016.dcf) | NBPSeq | 75601 | 80640 | | Unassigned |
| [017](config-017.dcf) | NBPSeq | 80641 | 85680 | | Unassigned |
| [018](config-018.dcf) | NBPSeq | 85681 | 90720 | | Unassigned |
| [019](config-019.dcf) | NBPSeq | 90721 | 95760 | | Unassigned |
| [020](config-020.dcf) | NBPSeq | 95761 | 100800 | | Unassigned |
| [021](config-021.dcf) | NBPSeq | 100801 | 105840 | | Unassigned |
| [022](config-022.dcf) | NBPSeq | 105841 | 110880 | | Unassigned |
| [023](config-023.dcf) | NBPSeq | 110881 | 115920 | | Unassigned |
| [024](config-024.dcf) | NBPSeq | 115921 | 120960 | | Unassigned |
| [025](config-025.dcf) | NBPSeq | 120961 | 126000 | | Unassigned |
| [026](config-026.dcf) | NBPSeq | 126001 | 131040 | | Unassigned |
| [027](config-027.dcf) | NBPSeq | 131041 | 136080 | | Unassigned |
| [028](config-028.dcf) | DESeq2 | 1 | 5040 | | Unassigned |
| [029](config-029.dcf) | DESeq2 | 5041 | 10080 | | Unassigned |
| [030](config-030.dcf) | DESeq2 | 10081 | 15120 | | Unassigned |
| [031](config-031.dcf) | DESeq2 | 15121 | 20160 | | Unassigned |
| [032](config-032.dcf) | DESeq2 | 20161 | 25200 | | Unassigned |
| [033](config-033.dcf) | DESeq2 | 25201 | 30240 | | Unassigned |
| [034](config-034.dcf) | DESeq2 | 30241 | 35280 | | Unassigned |
| [035](config-035.dcf) | DESeq2 | 35281 | 40320 | | Unassigned |
| [036](config-036.dcf) | DESeq2 | 40321 | 45360 | | Unassigned |
| [037](config-037.dcf) | DESeq2 | 45361 | 50400 | | Unassigned |
| [038](config-038.dcf) | DESeq2 | 50401 | 55440 | | Unassigned |
| [039](config-039.dcf) | DESeq2 | 55441 | 60480 | | Unassigned |
| [040](config-040.dcf) | DESeq2 | 60481 | 65520 | | Unassigned |
| [041](config-041.dcf) | DESeq2 | 65521 | 70560 | | Unassigned |
| [042](config-042.dcf) | DESeq2 | 70561 | 75600 | | Unassigned |
| [043](config-043.dcf) | DESeq2 | 75601 | 80640 | | Unassigned |
| [044](config-044.dcf) | DESeq2 | 80641 | 85680 | | Unassigned |
| [045](config-045.dcf) | DESeq2 | 85681 | 90720 | | Unassigned |
| [046](config-046.dcf) | DESeq2 | 90721 | 95760 | | Unassigned |
| [047](config-047.dcf) | DESeq2 | 95761 | 100800 | | Unassigned |
| [048](config-048.dcf) | DESeq2 | 100801 | 105840 | | Unassigned |
| [049](config-049.dcf) | DESeq2 | 105841 | 110880 | | Unassigned |
| [050](config-050.dcf) | DESeq2 | 110881 | 115920 | | Unassigned |
| [051](config-051.dcf) | DESeq2 | 115921 | 120960 | | Unassigned |
| [052](config-052.dcf) | DESeq2 | 120961 | 126000 | | Unassigned |
| [053](config-053.dcf) | DESeq2 | 126001 | 131040 | | Unassigned |
| [054](config-054.dcf) | DESeq2 | 131041 | 136080 | | Unassigned |
| [055](config-055.dcf) | DESeq | 1 | 13608 | | Unassigned |
| [056](config-056.dcf) | DESeq | 13609 | 27216 | | Unassigned |
| [057](config-057.dcf) | DESeq | 27217 | 40824 | | Unassigned |
| [058](config-058.dcf) | DESeq | 40825 | 54432 | | Unassigned |
| [059](config-059.dcf) | DESeq | 54433 | 68040 | | Unassigned |
| [060](config-060.dcf) | DESeq | 68041 | 81648 | | Unassigned |
| [061](config-061.dcf) | DESeq | 81649 | 95256 | | Unassigned |
| [062](config-062.dcf) | DESeq | 95257 | 108864 | | Unassigned |
| [063](config-063.dcf) | DESeq | 108865 | 122472 | | Unassigned |
| [064](config-064.dcf) | DESeq | 122473 | 136080 | | Unassigned |
| [065](config-065.dcf) | edgeR.lrt | 1 | 34020 | | Unassigned |
| [066](config-066.dcf) | edgeR.lrt | 34021 | 68040 | | Unassigned |
| [067](config-067.dcf) | edgeR.lrt | 68041 | 102060 | | Unassigned |
| [068](config-068.dcf) | edgeR.lrt | 102061 | 136080 | | Unassigned |
| [069](config-069.dcf) | edgeR.qlf | 1 | 34020 | | Unassigned |
| [070](config-070.dcf) | edgeR.qlf | 34021 | 68040 | | Unassigned |
| [071](config-071.dcf) | edgeR.qlf | 68041 | 102060 | | Unassigned |
| [072](config-072.dcf) | edgeR.qlf | 102061 | 136080 | | Unassigned |
| [073](config-073.dcf) | DSS | 1 | 34020 | | Unassigned |
| [074](config-074.dcf) | DSS | 34021 | 68040 | | Unassigned |
| [075](config-075.dcf) | DSS | 68041 | 102060 | | Unassigned |
| [076](config-076.dcf) | DSS | 102061 | 136080 | | Unassigned |
| [077](config-077.dcf) | voom | 1 | 34020 | | Unassigned |
| [078](config-078.dcf) | voom | 34021 | 68040 | | Unassigned |
| [079](config-079.dcf) | voom | 68041 | 102060 | | Unassigned |
| [080](config-080.dcf) | voom | 102061 | 136080 | | Unassigned |
| [081](config-081.dcf) | Wilcoxon | 1 | 34020 | | Unassigned |
| [082](config-082.dcf) | Wilcoxon | 34021 | 68040 | | Unassigned |
| [083](config-083.dcf) | Wilcoxon | 68041 | 102060 | | Unassigned |
| [084](config-084.dcf) | Wilcoxon | 102061 | 136080 | | Unassigned |
| [085](config-085.dcf) | T.test | 1 | 34020 | | Unassigned |
| [086](config-086.dcf) | T.test | 34021 | 68040 | | Unassigned |
| [087](config-087.dcf) | T.test | 68041 | 102060 | | Unassigned |
| [088](config-088.dcf) | T.test | 102061 | 136080 | | Unassigned |
| [089](config-089.dcf) | NOISeq | 1 | 34020 | | Unassigned |
| [090](config-090.dcf) | NOISeq | 34021 | 68040 | | Unassigned |
| [091](config-091.dcf) | NOISeq | 68041 | 102060 | | Unassigned |
| [092](config-092.dcf) | NOISeq | 102061 | 136080 | | Unassigned |
| [093](config-093.dcf) | ABSSeq | 1 | 34020 | | Unassigned |
| [094](config-094.dcf) | ABSSeq | 34021 | 68040 | | Unassigned |
| [095](config-095.dcf) | ABSSeq | 68041 | 102060 | | Unassigned |
| [096](config-096.dcf) | ABSSeq | 102061 | 136080 | | Unassigned |
| [097](config-097.dcf) | ROTS | 1 | 34020 | | Unassigned |
| [098](config-098.dcf) | ROTS | 34021 | 68040 | | Unassigned |
| [099](config-099.dcf) | ROTS | 68041 | 102060 | | Unassigned |
| [100](config-100.dcf) | ROTS | 102061 | 136080 | | Unassigned |
