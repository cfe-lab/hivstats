
This repository holds code used to analyze some statistical properties of HIV genome related to defects detection.

Results produced with it are used in CFEIntact.

# Methodology: Progressive Intactness Filtering

A key challenge in deriving statistical thresholds for HIV intactness is avoiding circular dependencies. We need to determine what size, distance, and indel impact values are "normal" for intact sequences, but we need to know which sequences are intact to calculate those distributions.

To solve this, we use **progressive intactness filtering** based on CFEIntact's own defect classifications:

## Level 1: Structural Intactness (for size analysis)

Used when deriving size distributions. Only checks CFEIntact defects that are independent of size/length:
- **APOBECHypermutation**: G-to-A hypermutation signature
- **NonHIV**: Sequence contains non-HIV DNA
- **Scramble**: Sequence order is scrambled
- **InternalInversion**: Part of sequence is inverted
- **UnknownNucleotide**: Contains ambiguous bases
- **MissingORF**: Expected ORF is missing
- **LongDeletion**: Large deletions in the sequence

These defects can be detected without knowing expected sizes or alignment distances.

## Level 2: Distance-Based Intactness (for distance analysis)

Used when deriving distance distributions. Includes Level 1 defects plus defects that depend on alignment quality:
- All structural defects from Level 1
- **Plus**:
  - **Deletion**: Deletion mutations
  - **Insertion**: Insertion mutations
  - **MutatedStartCodon**: Start codon is mutated
  - **MutatedStopCodon**: Stop codon is mutated
  - **InternalStop**: Internal stop codons

These defects relate to sequence alignment and codon integrity but don't depend on knowing expected distances.

## Level 3: Indel-Based Intactness (for indel impact analysis)

Used when deriving indel impact distributions. Includes structural defects plus actual indel defects:
- All defects from Level 1 (structural defects)
- **Plus**:
  - **Deletion**: Deletion mutations in the sequence
  - **Insertion**: Insertion mutations in the sequence

This focused approach only excludes sequences with clear structural problems or actual insertion/deletion defects, allowing us to analyze the impact of indels without being confounded by other distance-related metrics like sequence divergence, frameshifts, or stop codons (which may themselves be consequences of indels we're trying to measure).

## Why This Matters

This approach ensures we don't use metric-derived thresholds to define the very populations used to derive those thresholds. Each level builds on previous levels without creating circular dependencies.

**All intactness criteria come from CFEIntact itself** - we don't invent additional thresholds. We simply categorize CFEIntact's existing defect codes by what they depend on.

# How to use

First run

```shell
sh src/install-dependencies.sh
make all
```

Then run

```shell
make serve    # This will open Jupyter notebook with the results.
```

# Inputs

Files in `inputs` directory contain all data used in the analysis.

Below is a description of individual files:

## `los-alamos-all-sequences.fasta`

All sequences of subtype B in the Los-Alamos database.

## `los-alamos-plasma-sequences.fasta`

Sequences of subtype B that were extracted from plasma.

Downloaded from Los-Alamos database.

## `other-intact-sequences.fasta`

Subtype B sequences from CFEIntact's database.

Not really required, but nice to have since they are definitely intact, but weren't in the other `.fasta` files.

## `individual-plasma/seq/*`

Also downloaded from the Los-Alamos database, extracted from plasma.

- These are subtype B
- 1 per patient
- Clipped to the region
