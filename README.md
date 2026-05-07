# NgoAMR: Identification of genetic determinants of resistance in *Neisseria gonorrhoeae* using ARIBA

A repository of tools and databases for the identification of genetic determinants of antimicrobial resistance (AMR) in the human pathogen *Neisseria gonorrhoeae* (Ngo) from Illumina short-read whole-genome sequencing (WGS) data.

The core tool for resistance identification is [ARIBA](https://sanger-pathogens.github.io/ariba/) (Antibiotic Resistance Identification By Assembly), which performs reference-guided local assembly of target loci directly from read pairs. This repository provides a batch-processing wrapper script and a custom *N. gonorrhoeae* ARIBA reference database.

---

## Repository structure

```
NgoAMR/
├── tools/
│   └── ngoamr_ariba_batch_v0.3.py     # Batch ARIBA runner and report processor
├── databases/
│   └── ARIBA_NgoDB_070526/             # Compiled ARIBA reference database
└── data/
    └── input_files_ARIBAdb/            # Source sequences and metadata used to build the database
```

---

## Requirements

- [ARIBA](https://sanger-pathogens.github.io/ariba/) v2.14.6 or compatible, with its dependencies:
  - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) ≥ 2.3.5
  - [CD-HIT](https://sites.google.com/view/cd-hit) ≥ 4.8
  - [MUMmer](http://mummer.sourceforge.net/) ≥ 3.1
- Python ≥ 3.7

---

## Tools

### `ngoamr_ariba_batch_v0.3.py`

A Python script that automates running ARIBA over multiple samples spread across one or more directories, post-processes the individual ARIBA reports to enable correct detection of two specific indels, and builds a combined summary table across all samples.

**Workflow:**
1. Scans the specified input directories for paired-end `.fastq.gz` files 
2. For each sample, runs `ariba run` against the specified database
3. Patches the `report.tsv` output to produce a `report_complete.tsv` that reports the **penA insD345 insertion** and the **mtrR promoter −53A deletion** (see [Special indel handling](#special-indel-handling) below)
4. Runs `ariba summary` over all `report_complete.tsv` files to produce a combined resistance matrix

**Supported FASTQ naming conventions:** `_1`/`_2.fastq.gz`, `_R1`/`_R2.fastq.gz`, `_R1_001`/`_R2_001.fastq.gz`, `_forward`/`_reverse.fq.gz`

**Smart re-run logic:** samples with an existing `report_complete.tsv` are skipped automatically. Incomplete output folders (missing `report.tsv`) are deleted and re-run from scratch.

#### Usage

```bash
python ngoamr_ariba_batch_v0.3.py -d <INPUT_DIR(S)> -p <DB_PATH> -o <OUTPUT_DIR> [options]
```

#### Arguments

| Flag | Long form | Required | Default | Description |
|------|-----------|:--------:|:-------:|-------------|
| `-d` | `--input_dirs` | Yes | — | Comma-separated list of directories containing `.fastq.gz` files |
| `-p` | `--db_path` | Yes | — | Path to the compiled ARIBA database (e.g. `databases/ARIBA_NgoDB_070526`) |
| `-o` | `--output_dir` | Yes | — | Output directory for per-sample ARIBA results and the combined summary |
| `-t` | `--threads` | No | 1 | Number of threads passed to `ariba run` |
| `-m` | `--modify_report` | No | False | Only generate `report_complete.tsv` from an existing `report.tsv` (skip `ariba run`). Useful for recovering from partial runs or standalone `ariba` results|
| `-s` | `--summary_only` | No | False | Skip `ariba run` entirely; only re-run `ariba summary` over already-processed samples |
| `-n` | `--short_names` | No | False | Use folder names as sample labels in the summary instead of full paths. Note: full paths are required for compatibility with downstream tools such as sensityping |

#### Example

```bash
python ngoamr_ariba_batch_v0.3.py \
  -d /data/run1_fastqs,/data/run2_fastqs \
  -p databases/ARIBA_NgoDB_070526 \
  -o results/ariba_output \
  -t 4
```

#### Special indel handling

ARIBA does not natively propagate certain indels through to the `ariba summary` output. The script automatically patches two fields in each `report.tsv` before building the summary to enable detection of:

- **penA insD345 insertion** — identified by the presence of `D147_T148insT` or `R146_D147insR` in the report
- **mtrR promoter −53A deletion** — identified by the presence of `A197.` in the report

The patched file is saved as `report_complete.tsv` alongside the original `report.tsv`.

---

## Databases

### ARIBA_NgoDB_070526

A custom ARIBA reference database for *N. gonorrhoeae* AMR detection, compiled from 1,082 curated sequences covering 35 gene groups (1,075 coding sequences and 7 non-coding promoter sequences across 3 promoter loci). Built with `ariba prepareref` as follows:

```bash
ariba prepareref -f sequences_070526.fa -m metadata_070526.tsv --cdhit_clusters clusters_070526.tsv ARIBA_NgoDB_040526
```

Source files are in `data/input_files_ARIBAdb/`.
- `sequences_070526.fa` includes all reference sequences, including multiple versions for 8 gene groups to represent diverse alleles and recombinant variants: mtrD (118 alleles, including 3 mosaic sequences), penA (921 alleles), mtrR_promoter (4 alleles, including 3 mosaic sequences), porB1b (4 alleles), and folP, norM_promoter, ponA, and porA (2 alleles each)
- `metadata_070526.tsv` includes information on known mutations or gene presence associated with AMR for all the sequences in the database. Each row defines one entry; sequences with multiple tracked variants have one row per variant. The six tab-separated columns are:

  | Column | Field | Description |
  |--------|-------|-------------|
  | 1 | Sequence name | Must match the FASTA header exactly |
  | 2 | Gene status | `1` = coding sequence; `0` = non-coding sequence |
  | 3 | Sequence type | `0` = presence/absence detection (no variant calling); `1` = variant calling enabled |
  | 4 | Variant | Amino acid change for coding sequences (e.g. `R228S`) or nucleotide change for non-coding (e.g. `C1184T`); `.` if none |
  | 5 | Variant group | Label used to group related variants across sequences in the summary output; `.` if unused |
  | 6 | Description | Free-text annotation; `.` if unused |
- `clusters_070526.tsv` predefines a single cluster for each of the 8 multi-allele gene groups listed above, to prevent CD-HIT from splitting recombinant or divergent alleles into separate clusters


#### Genes and mutations

The table below lists all genes and mutations tracked in the database, organised by antimicrobial class. For genes with many reference alleles (*penA*, *mtrD*), key variant positions are summarised rather than listed exhaustively.

**Detection types:**
- *Allele* — gene detected by best-matching reference allele; variants called against that allele
- *Presence/absence* — gene flagged as present or absent; no variant calling

---

##### Extended-spectrum cephalosporins (ESCs) and beta-lactams

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *penA* | Penicillin-binding protein 2 (PBP2) | Variant | Variants at positions A311, T316, T483S, A501V/P/T, N512Y, T534A, G542S, G545S, P551S/L, insD345 insertion | insD345 requires special indel handling (`penA100partial` anchor sequence used for identification) |
| *ponA* | Penicillin-binding protein 1 (PBP1) | Variant | L421P | High-level chromosomal penicillin resistance |
| *rpoB* | RNA polymerase β subunit | Variant | P157L, G158V, R201H, H553N | P157L/G158V/R201H reduce ESC susceptibility; H553N confers rifampicin resistance |
| *rpoD* | RNA polymerase σ factor | Variant | E98K | Decreased susceptibility to ESCs |
| *blaTEM* | TEM β-lactamase | Presence/absence | — | Acquired; confers penicillin resistance |

---

##### Macrolides (azithromycin)

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *23S* rDNA | 23S ribosomal RNA | Variant | A2045G, C2597T | A2045G (E. coli A2059G) high-level resistance; C2597T (E. coli C2611T) low-level resistance |
| *mef* | MefE efflux pump | Presence/absence | — | Acquired macrolide efflux |
| *ereA* | EreA macrolide esterase | Presence/absence | — | Acquired macrolide inactivation |
| *ereB* | EreB macrolide esterase | Presence/absence | — | Acquired macrolide inactivation |
| *ermA* | ErmA rRNA methylase | Presence/absence | — | Acquired; methylates 23S rRNA |
| *ermB* | ErmB rRNA methylase | Presence/absence | — | Acquired; methylates 23S rRNA |
| *ermC* | ErmC rRNA methylase | Presence/absence | — | Acquired; methylates 23S rRNA |
| *ermF* | ErmF rRNA methylase | Presence/absence | — | Acquired; methylates 23S rRNA |
| *mtrD* | MtrCDE efflux pump inner membrane protein | Variant and allele typing | R714G, K823E | All alleles tracked to identify mosaicism; mutations increase azithromycin MIC |
| *mtrR* | MtrR transcriptional repressor | Variant | A39T, G45D | Loss-of-function upregulates MtrCDE efflux |
| *mtrR*_promoter | *mtrR* promoter region | Variant | G120A (mtr-120 mutation), A195C (−38A>C); mosaic alleles 1–3 | Promoter mutations upregulate MtrCDE efflux; −53A deletion requires special indel handling (A197.) |
| *macAB*_promoter | *macAB* promoter region | Variant | G178T (−48G>T) | Overexpression of MacAB efflux pump |
| *rplD* | 50S ribosomal protein L4 | Variant | G70D | Decreased ribosome binding; increased azithromycin MIC |
| *mtrC* | MtrCDE outer membrane lipoprotein | Presence/absence | — | Loss-of-function may affect azithromycin susceptibility |

---

##### Fluoroquinolones

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *gyrA* | DNA gyrase subunit A (GyrA) | Variant | S91F, S91I, S91T, S91Y, A92T, D95G, D95N | S91 and D95 mutations confer fluoroquinolone resistance; A92T associated with gepotidacin resistance |
| *gyrB* | DNA gyrase subunit B (GyrB) | Variant | Any D429 or K450 substitutions | All possible amino acid substitutions at each position; associated with zoliflodacin resistance |
| *parC* | Topoisomerase IV subunit C (ParC) | Variant | D86N, S87A, S87I, S87N, S87R, S87W, S88P, E91K | Contribute to fluoroquinolone resistance |
| *parE* | Topoisomerase IV subunit E (ParE) | Variant | G410V | Contributes to fluoroquinolone resistance |
| *norM*_promoter | *norM* promoter region | Variant | A328G (−7A>G, RBS mutation), C231T (−104C>T, −35 box) | Both mutations cause overexpression of the NorM efflux pump |

---

##### Tetracyclines

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *tetM* | Tet(M) ribosomal protection protein | Presence/absence | — | Acquired; plasmid-borne tetracycline resistance |
| *rpsJ* | 30S ribosomal protein S10 | Variant | V57L, V57M, V57Q | Chromosomal tetracycline resistance |

---

##### Spectinomycin

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *16S* rDNA | 16S ribosomal RNA | Variant | C1184T | E. coli C1192T equivalent; confers spectinomycin resistance |
| *rpsE* | 30S ribosomal protein S5 | Variant | T24P, K28E | Disrupt spectinomycin binding to ribosomal target |

---

##### Sulfonamides

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *folP* | Dihydropteroate synthase (FolP) | Variant | R228S | Confer sulfonamide resistance |

---

##### Outer membrane permeability (porins)

| Gene | Gene product | Detection | Key mutations / variants tracked | Notes |
|------|-------------|-----------|----------------------------------|-------|
| *porA* | Outer membrane porin PorA | Presence/absence | — | Used for porin typing |
| *porB1b* | Outer membrane porin PorB1b | Variant | G120K, A121D/N | Reduce porin permeability; contribute to decreased antimicrobial uptake |
| *porB1a* | Outer membrane porin PorB1a | Presence/absence | — | Used for porin typing |
| *porB2* | Outer membrane porin PorB2 | Presence/absence | — | Used for porin typing |
| *porB3* | Outer membrane porin PorB3 | Presence/absence | — | Used for porin typing |

---

## Output

For each sample, the script creates a subdirectory `<sample_name>_ARIBA/` inside the output directory, containing the standard ARIBA output files plus:

- `report_complete.tsv` — modified version of `report.tsv` with patched indel fields enabling correct detection of *penA* insD345 and *mtrR* promoter −53A deletion in the summary step

After processing all samples, the following summary files are written to the working directory:

- `ariba_summary.csv` — binary presence/absence and variant matrix across all samples
- `ariba_summary.phandango.csv` / `.tre` — files for visualisation in [Phandango](https://jameshadfield.github.io/phandango/)
- `filenames.txt` — list of `report_complete.tsv` paths used as input to `ariba summary`

---

## Author

Leonor Sánchez Busó (last updated 06/05/2026)
