
This repository holds code used to analyze some statistical properties of HIV genome related to defects detection.

Results produced with it are used in CFEIntact.

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

---

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

---

# Kernel Density Estimation (KDE)

All distribution plots now include **Kernel Density Estimation** curves overlaid on the histograms. KDE provides a smooth, continuous estimate of the underlying probability density without assuming any particular distribution shape (like Normal or unimodal).

## How KDE Works

KDE places a "bump" (Gaussian kernel) on each data point and sums all bumps together:

$$\hat{f}(x) = \frac{1}{n h}\sum_{i=1}^n K\left(\frac{x-x_i}{h}\right)$$

where:
- $K(u)$ is the Gaussian kernel: $\frac{1}{\sqrt{2\pi}}e^{-u^2/2}$
- $h$ is the bandwidth (controls smoothness)
- $x_i$ are the observed data points

## Bandwidth Selection

We use **Silverman's rule of thumb** to automatically select bandwidth:

$$h \approx 0.9 \cdot \min(\sigma, \text{IQR}/1.34) \cdot n^{-1/5}$$

This balances smoothness and detail:
- **Small bandwidth**: Shows more detail, may reveal multiple peaks
- **Large bandwidth**: Smoother curve, may merge nearby peaks

## Why KDE is Useful

1. **No distributional assumptions**: Unlike parametric methods, KDE doesn't assume data is Normal
2. **Reveals multiple peaks**: Can show bimodal or multimodal distributions naturally
3. **Visual validation**: Helps identify whether outlier filtering is appropriate
4. **Continuous estimates**: Provides smooth density for any metric value

The KDE curves appear as red lines on single-distribution plots, and as dashed lines (dark blue for intact, dark red for defective) on comparison plots.

---
