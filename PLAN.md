# Plan: Align 13 MRSA assemblies against USA300 FPR3757 with minimap2

## Context

The user has downloaded 13 gzipped *S. aureus* assemblies into `assemblies/assembly1.fasta.gz`–`assembly13.fasta.gz` (each ~850 KB compressed, ~3 contigs per assembly). They need to align each assembly against the "official MRSA reference at NCBI" and have selected **USA300 FPR3757 — RefSeq `NC_007793.1`** (community-acquired MRSA, complete genome). The outputs will be sorted/indexed BAM files plus per-sample coverage stats, which will feed downstream variant/SNP analysis.

**Notable discrepancy:** the repo's `README.md` currently says "MRSA isolate alignment against NCTC 8325 reference". Since the user chose USA300 FPR3757, the README should be updated to match (or left to the user to update separately).

**Environment findings from exploration:**
- `minimap2` and `samtools` are **NOT installed** — must be installed first.
- `mamba` (2.5.0) + `conda` (26.1.1) are available via miniforge3; currently in `base` env.
- `efetch` (entrez-direct 24.7) is available at `/home/biouser/edirect/` for the NCBI download.
- 8 CPUs, ~937 GB free disk — ample headroom.
- No existing `reference/` dir and no existing scripts/pipelines in the repo.

## Approach

### 1. Create a dedicated conda env for alignment tools

Keep `base` clean; install minimap2 + samtools into a new env named `mrsa-align`:

```bash
mamba create -y -n mrsa-align -c bioconda -c conda-forge minimap2 samtools
mamba activate mrsa-align
```

Expected versions: minimap2 ≥ 2.26, samtools ≥ 1.19.

### 2. Fetch the USA300 FPR3757 reference

```bash
mkdir -p reference alignments logs
efetch -db nuccore -id NC_007793.1 -format fasta > reference/USA300_FPR3757.fasta
samtools faidx reference/USA300_FPR3757.fasta
```

Sanity check: the `.fai` should show one sequence of ~2,872,769 bp (chromosome). Note: `NC_007793.1` is the chromosome only; the USA300 FPR3757 plasmids (`NC_007790–NC_007792`) are excluded unless the user wants them added.

### 3. Align each assembly with minimap2 `-x asm5`

These are same-species assemblies (~99%+ ANI expected against USA300), so the `asm5` preset is correct. Loop over all 13 assemblies, pipe directly into `samtools sort`, then index:

```bash
REF=reference/USA300_FPR3757.fasta
for i in $(seq 1 13); do
  minimap2 -ax asm5 -t 6 "$REF" "assemblies/assembly${i}.fasta.gz" \
    2> "logs/assembly${i}.mm2.log" \
  | samtools sort -@ 2 -o "alignments/assembly${i}.bam" -
  samtools index "alignments/assembly${i}.bam"
done
```

Notes:
- `-a` emits SAM (required for downstream BAM).
- `-t 6` for minimap2 + `-@ 2` for sort ≈ full use of 8 cores without oversubscribing.
- minimap2 reads gzipped FASTA natively — no need to decompress.

### 4. Per-sample QC summary

Collect mapping stats into a single TSV so outliers are easy to spot:

```bash
{
  printf "sample\tmapped_contigs\ttotal_contigs\tmean_depth\tbreadth_pct\n"
  for i in $(seq 1 13); do
    bam="alignments/assembly${i}.bam"
    mapped=$(samtools view -c -F 260 "$bam")   # primary mapped
    total=$(samtools view -c -F 256 "$bam")    # primary records
    read depth breadth < <(samtools coverage "$bam" \
      | awk 'NR>1 {d+=$7*$3; l+=$3; c+=$6*$3} END {printf "%.2f %.2f", d/l, c/l}')
    printf "assembly%d\t%d\t%d\t%s\t%s\n" "$i" "$mapped" "$total" "$depth" "$breadth"
  done
} > alignments/coverage_summary.tsv
```

Flag any assembly with `breadth_pct < 95` or an unusually high fraction of unmapped contigs — those are candidates for a rerun with `-x asm10` (more divergence tolerance) or closer inspection.

### 5. Save the run as a script (optional but recommended)

Persist the pipeline in `scripts/align_assemblies.sh` so the run is reproducible. Exact shell body above; add `set -euo pipefail` and an activation check for the `mrsa-align` env.

### 6. Commit `PLAN.md` and push to GitHub

Before running the pipeline, save this plan as `PLAN.md` at the repo root so it's tracked in version control, then commit and push to `origin/main` (`github.com/Bri-W12/MRSA-analysis`):

```bash
cp /home/biouser/.claude/plans/i-need-to-align-misty-kettle.md /home/biouser/MRSA-analysis/PLAN.md
cd /home/biouser/MRSA-analysis
git add PLAN.md
git commit -m "Add alignment plan for USA300 FPR3757 reference"
git push origin main
```

Note: the `assemblies/` directory (~11 MB of gzipped FASTAs) is currently untracked. It will **not** be included in the PLAN commit — only `PLAN.md` is staged explicitly. If the assemblies, BAMs, or reference should also be version-controlled (or excluded via `.gitignore`), flag this separately.

## Files to be created

| Path | Purpose | Committed to git? |
|------|---------|---|
| `PLAN.md` | This plan, at repo root | **Yes** |
| `scripts/align_assemblies.sh` | Reproducible pipeline script | Not in this commit |
| `reference/USA300_FPR3757.fasta` (+ `.fai`) | Reference genome + index | No (regenerable) |
| `alignments/assembly{1..13}.bam` (+ `.bai`) | Sorted, indexed BAMs | No (regenerable) |
| `alignments/coverage_summary.tsv` | Per-sample QC table | No (regenerable) |
| `logs/assembly{1..13}.mm2.log` | minimap2 stderr per sample | No |

No existing project files are modified. (The `README.md` currently references NCTC 8325; flagged for the user to update separately.)

## Verification

1. **Environment:** `mamba activate mrsa-align && minimap2 --version && samtools --version | head -1` — both should print versions without errors.
2. **Reference integrity:** `samtools faidx reference/USA300_FPR3757.fasta && cat reference/USA300_FPR3757.fasta.fai` — expect a single ~2.87 Mb contig.
3. **Run the pipeline:** `bash scripts/align_assemblies.sh` — completes in a few minutes on 8 cores for 13 × ~3-contig assemblies.
4. **Per-BAM smoke test:** `samtools quickcheck alignments/*.bam && echo OK` — all BAMs must pass.
5. **QC table review:** `column -t alignments/coverage_summary.tsv` — all 13 samples should show breadth ≳ 95% and mean depth ≈ 1 (since we're aligning one set of contigs per sample). Any sample materially below 95% breadth should be rerun with `-x asm10` to check whether it's truly divergent or just preset-sensitive.
6. **Spot-check visualization (optional):** load one BAM + reference into IGV to confirm contigs tile the reference end-to-end.
7. **Git push confirmation:** after `git push origin main`, verify the commit appears at `https://github.com/Bri-W12/MRSA-analysis/blob/main/PLAN.md`.
