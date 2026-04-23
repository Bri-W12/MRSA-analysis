#!/usr/bin/env bash
set -euo pipefail

# Call SNPs/indels for 13 MRSA assemblies vs USA300 FPR3757 using bcftools mpileup
# Requires: mrsa-align conda env (minimap2 2.30, samtools 1.23.1, bcftools 1.23.1)
# Run align_assemblies.sh first to generate BAM files
# Usage: micromamba run -n mrsa-align bash scripts/call_variants.sh

REF=reference/USA300_FPR3757.fasta
mkdir -p variants

for i in $(seq 1 13); do
    bam="alignments/assembly${i}.bam"

    bcftools mpileup -Ou -f "$REF" "$bam" \
        --min-MQ 30 --min-BQ 20 \
    | bcftools call -mv -Oz -o "variants/assembly${i}.vcf.gz"
    bcftools index "variants/assembly${i}.vcf.gz"
done

# Merge into a single multi-sample VCF
bcftools merge -Oz -o variants/all_assemblies.vcf.gz \
    variants/assembly{1..13}.vcf.gz
bcftools index variants/all_assemblies.vcf.gz

# Summary table
{
    printf "sample\tsnps\tindels\ttotal_variants\n"
    for i in $(seq 1 13); do
        vcf="variants/assembly${i}.vcf.gz"
        snps=$(bcftools view -v snps "$vcf" 2>/dev/null | grep -vc "^#")
        indels=$(bcftools view -v indels "$vcf" 2>/dev/null | grep -vc "^#")
        printf "assembly%d\t%d\t%d\t%d\n" "$i" "$snps" "$indels" "$((snps+indels))"
    done
} > variants/variant_summary.tsv

echo "Done. Results in variants/"
echo "  Per-sample VCFs:     variants/assembly{1..13}.vcf.gz"
echo "  Merged VCF:          variants/all_assemblies.vcf.gz"
echo "  Summary table:       variants/variant_summary.tsv"
