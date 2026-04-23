#!/usr/bin/env bash
set -euo pipefail

# Align 13 MRSA assemblies against USA300 FPR3757 (NC_007793.1) with minimap2
# Requires: mrsa-align conda env (minimap2 2.30, samtools 1.23.1)
# Usage: micromamba run -n mrsa-align bash scripts/align_assemblies.sh

REF=reference/USA300_FPR3757.fasta

mkdir -p reference alignments logs

# Fetch and index reference if not present
if [[ ! -f "$REF" ]]; then
    ~/edirect/efetch -db nuccore -id NC_007793.1 -format fasta > "$REF"
    samtools faidx "$REF"
fi

# Align each assembly; asm10 used for all (more permissive, better breadth for divergent isolates)
# Assemblies 11 and 12 are close USA300 relatives — asm5 would also work but asm10 is consistent
for i in $(seq 1 13); do
    minimap2 -ax asm10 -t 6 "$REF" "assemblies/assembly${i}.fasta.gz" \
        2> "logs/assembly${i}.mm2.log" \
    | samtools sort -@ 2 -o "alignments/assembly${i}.bam" -
    samtools index "alignments/assembly${i}.bam"
done

# QC summary
{
    printf "sample\tmapped_contigs\ttotal_contigs\tmean_depth\tbreadth_pct\n"
    for i in $(seq 1 13); do
        bam="alignments/assembly${i}.bam"
        mapped=$(samtools view -c -F 260 "$bam")
        total=$(samtools view -c -F 256 "$bam")
        read depth breadth < <(samtools coverage "$bam" \
            | awk 'NR>1 {d+=$7*$3; l+=$3; c+=$6*$3} END {printf "%.2f %.2f", d/l, c/l}')
        printf "assembly%d\t%d\t%d\t%s\t%s\n" "$i" "$mapped" "$total" "$depth" "$breadth"
    done
} > alignments/coverage_summary.tsv

echo "Done. QC summary: alignments/coverage_summary.tsv"
