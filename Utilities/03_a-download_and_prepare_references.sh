#!/bin/bash
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --mincpus 1
#SBATCH -o prepare_references.%j.out
#SBATCH -e prepare_references.%j.err
#SBATCH -J prepare_references

set -exuo pipefail

mkdir -p data/references
cd data/references

## a) Retrieve GRCh38.p13 references

fasta_url="http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_lcl="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz"
gtf_lcl="gencode.v41.primary_assembly.annotation.gtf.gz"

wget -O ${fasta_lcl} ${fasta_url} 
wget -O ${gtf_lcl} ${gtf_url}

## b) Modify the genome reference FASTA contig names
# Gencode use "chr1" etc versus Ensembl "1" etc.  Unplaced and unlocalized
# sequences such as "GL456210.1" have the same names in both versions.
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome

fasta_gencoded="Homo_sapiens.GRCh38.dna.primary_assembly.gencoded.fa"
zcat ${fasta_lcl} | 
    sed -E -e 's/^>(\S+).*/>\1 \1/' \
           -e 's/^>([0-9]+|[XY]) />chr\1 /' \
           -e 's/^>MT />chrM /' \
    > ${fasta_gencoded}

## c) Modify gene/treanscript/exon IDs in GTF
# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages.

ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
gtf_modified="gencode.v41.primary_assembly.annotation.mod_id.gtf.gz"
zcat ${gtf_lcl} |
    sed -E -e 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
           -e 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
           -e 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' |
    gzip -c - \
    > ${gtf_modified}

## d) Filter GTF

BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
zcat ${gtf_modified} \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > tmp_gene_allowlist

# Filter the GTF file based on the gene allowlist
gtf_filtered="gencode.v41.primary_assembly.annotation.filtered.gtf"
# Copy header lines beginning with "#"
zgrep -E "^#" ${gtf_modified} > ${gtf_filtered}
# Filter to the gene allowlist
zgrep -Ff tmp_gene_allowlist ${gtf_modified} >> ${gtf_filtered}

rm -f tmp_gene_allowlist
rm -f ${gtf_modified}
rm -f ${fasta_lcl}
rm -f ${gtf_lcl}
