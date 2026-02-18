# Consensus CLI (BCR/TCR)

CLI to compute consensus per `gene + exon` from IMGT/GENE-DB FASTA input.

## Features
- Parses IMGT FASTA headers (`accession|gene*allele|...|EXx|...`)
- Handles RTF-wrapped FASTA exports (filters control lines and trailing `\`)
- Groups sequences by base gene (`TRGC2` from `TRGC2*01`) and exon (`EX1`, `EX2`, `EX3`, etc.)
- Builds a reference-anchored multiple alignment for each gene+exon group
- Emits consensus with requested encoding:
  - `x`: unshared 5' or 3' ends
  - `X`: variable position (single base change or internal indel)
  - `acgtn`: shared base
- Writes one XLSX tab per `gene_exon` plus a summary tab
- Writes consensus FASTA entries per `gene|exon`
- Supports curated compiled CSV input with columns:
  - `gene`
  - `coordinate` or `coordinates` (e.g. `start..end+start..end`)
  - `region` (e.g. `EX1+EX2+EX3`)
  - `fasta`
  - if coordinates are missing for multi-exon rows, the tool can infer split lengths from `TCRC_exons.fasta` in the same folder (or `--exon-reference-fasta`)
- For compiled CSV runs, writes a `junctions` worksheet with:
  - split source (`coordinate_split`, `inferred_from_exons`, or `single_region_unsplit`)
  - inferred exon chain
  - 1-based junction boundaries in the combined sequence
- Writes a primary `gene_summary` worksheet with full-length per-gene consensus (exon-chain concatenation).
- Writes explicit isotype comparison worksheets per configured pair (`--compare-pair`).
  - defaults: `TRGC1,TRGC2` and `TRBC1,TRBC2`
  - includes mirrored row and pair consensus:
    - `N` for explicit base mismatches
    - `X` when ambiguity (`X`) is present in either input consensus
- Family collapse in `gene_summary` is strict by default for `EX2_family=EX2,EX2R,EX2T`:
  - collapse happens only when a gene has more than one of those members present.
  - custom families can be added/overridden via repeatable `--family FAMILY=MEM1,MEM2,...`.

## Install
```bash
cd /Users/troynoordenbos/code
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run
```bash
python consensus_cli.py \
  --input '/Users/troynoordenbos/code/xenium_fasta_BCR_TCR/input_data/TCRC_exons.fasta' \
  --output '/Users/troynoordenbos/code/xenium_fasta_BCR_TCR/consensus_regions.xlsx'
```

Compiled CSV input:
```bash
python consensus_cli.py \
  --input '/Users/troynoordenbos/code/xenium_fasta_BCR_TCR/input_data/TCR_compiled_fasta.csv' \
  --output '/Users/troynoordenbos/code/xenium_fasta_BCR_TCR/consensus_regions_compiled.xlsx'
```

Optional: restrict to specific genes/exons
```bash
python consensus_cli.py \
  --input '/path/to/TCRC_exons.fasta' \
  --output '/path/to/consensus_regions.xlsx' \
  --gene TRGC1 --gene TRGC2 \
  --exon EX2 --exon EX3 \
  --consensus-fasta-output '/path/to/consensus_regions.fasta'
```
