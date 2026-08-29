# PooledPP

**PooledPP: Pseudo-genome-based Pangenome Construction from Pooled HiFi Sequencing Data**

PooledPP is a computational pipeline developed for constructing a **pseudo-genome-based graph pangenome from pooled PacBio HiFi sequencing data**.

Unlike conventional pangenome approaches that rely on independently assembled genomes from multiple individuals, PooledPP extracts structural variant (SV) signatures directly from pooled HiFi read alignments, clusters redundant read-supported SV events, reconstructs representative alternative sequences, distributes compatible SVs among multiple pseudo-genomes, and subsequently incorporates these pseudo-genomes into a graph pangenome.

The resulting pseudo-genomes are computational representations of alternative sequence paths recovered from pooled long-read data and should **not** be interpreted as reconstructed genomes or phased haplotypes of individual animals.

> **Important:** The final graph-to-VCF conversion step requires **vg v1.40.0 or an earlier compatible version**.  
> **vg versions newer than v1.40.0 are not supported by the current PooledPP implementation.**

---

## Overview

Pooled HiFi sequencing provides long and highly accurate reads from multiple individuals, but individual-level genome identity is not retained. Consequently, conventional assembly-based pangenome construction is difficult to apply directly.

PooledPP addresses this problem by converting read-supported structural variation into explicit alternative sequence paths and organizing these paths into a set of pseudo-genomes for graph construction.

```text
Pooled HiFi reads
        |
        v
Alignment to reference genome
        |
        v
01.extractSignature.pl
        |
        |  SV signatures
        |  read-supported alternative sequences
        v
02.cluster.pl
        |
        |  clustered representative SV events
        v
03.generatePseudogenome.pl
        |
        |  Pseudogenome1.fa
        |  Pseudogenome2.fa
        |  ...
        v
Minigraph
        |
        |  graph pangenome
        v
GFA
        |
        v
04.generateVCF.pl
        |
        |  vg view
        |  vg deconstruct
        |  VCF normalization
        v
Graph-derived VCF
```

An optional BAM downsampling script is also provided to standardize sequencing depth when required.

---

## Key features

- Designed for **pooled PacBio HiFi sequencing data**
- Does not require individual-level de novo genome assemblies
- Extracts SV signatures directly from long-read alignments
- Uses breakpoint, split-alignment, CIGAR, and local anchor information
- Evaluates local sequencing depth around candidate SV anchors
- Supports simple and complex structural variation
- Clusters redundant SV signatures supported by multiple reads
- Retains representative read-derived alternative sequences
- Converts pooled SV evidence into sequence-resolved pseudo-genomes
- Separates nearby or conflicting SV events among different pseudo-genomes
- Constructs a graph pangenome using Minigraph
- Deconstructs the graph into VCF representation using vg
- Provides an optional BAM downsampling utility

---

## Repository structure

```text
PooledPP/
├── 00.downSamplingBam.sh
├── 01.extractSignature.pl
├── 02.cluster.pl
├── 03.generatePseudogenome.pl
├── 04.generateVCF.pl
└── README.md
```

---

## Requirements

### External software

PooledPP requires the following programs:

```text
Perl
samtools
bedtools
csvtk
minigraph
vg
bcftools
bgzip
tabix
gzip
```

The optional BAM downsampling script additionally requires:

```text
awk
bc
```

### Perl modules

The following Perl modules are used:

```text
Getopt::Long
IO::Compress::Gzip
Set::IntervalTree
File::Temp
Data::Dumper
```

Missing Perl modules can be installed using CPAN or Conda.

For example:

```bash
cpan Set::IntervalTree
```

---

## IMPORTANT: vg version requirement

> **PooledPP requires vg v1.40.0 or an earlier compatible release for the graph-to-VCF conversion step.**

The current implementation of `04.generateVCF.pl` relies on the vg workflow used in **vg v1.40.0 and earlier compatible versions**.

The script directly performs:

```bash
vg view -F -g input.gfa > output.vg
vg deconstruct output.vg -a > output.vcf
```

For reproducibility, we strongly recommend:

```text
vg v1.40.0
```

The supported workflow is:

```text
vg <= v1.40.0
```

Versions newer than vg v1.40.0 are **not supported by the current PooledPP implementation**, because the graph-to-VCF conversion procedure required by `04.generateVCF.pl` is no longer available in the same form in later vg releases.

Using a newer vg version may therefore cause the final graph-to-VCF conversion step to fail.

Before running `04.generateVCF.pl`, check the installed vg version:

```bash
vg version
```

or, depending on the installation:

```bash
vg --version
```

The recommended version is:

```text
vg v1.40.0
```

---

## Input requirements

The main PooledPP workflow requires:

1. A reference genome in FASTA format
2. Pooled HiFi reads aligned to the reference genome
3. A coordinate-sorted BAM file
4. A read-name-sorted BAM file

The reference genome should be indexed using:

```bash
samtools faidx reference.fa
```

The BAM alignment should retain the alignment and sequence information required for reconstructing read-supported SV paths and alternative sequences.

---

# Workflow

## Step 0. Optional BAM downsampling

`00.downSamplingBam.sh` can be used to normalize sequencing depth before SV signature extraction.

### Usage

```bash
bash 00.downSamplingBam.sh \
    input.bam \
    output.bam \
    target_depth
```

### Example

```bash
bash 00.downSamplingBam.sh \
    sample.bam \
    sample.10x.bam \
    10
```

The script:

1. calculates the average sequencing depth of the input BAM;
2. calculates the required sampling fraction;
3. downsamples the BAM using `samtools view`;
4. indexes the resulting BAM; and
5. reports the sequencing depth after downsampling.

Outputs:

```text
sample.10x.bam
sample.10x.bam.bai
```

This step is optional.

---

## Step 1. Extract structural variant signatures

SV signatures are extracted using:

```text
01.extractSignature.pl
```

Two representations of the same BAM file are required:

- coordinate-sorted BAM;
- read-name-sorted BAM.

### Prepare BAM files

```bash
samtools sort \
    -o sample.coordinate.sorted.bam \
    sample.bam

samtools sort -n \
    -o sample.readname.sorted.bam \
    sample.bam

samtools index sample.coordinate.sorted.bam
```

### Run signature extraction

```bash
perl 01.extractSignature.pl \
    --bam sample.coordinate.sorted.bam \
    --bamSortRead sample.readname.sorted.bam \
    --output sample.signature
```

The output is automatically compressed:

```text
sample.signature.gz
```

### Major default parameters

```text
--min-mapq                 20
--min-align-len            1000
--min-sv-len               50
--max-sv-len               10000

--anchor-len-insertion     230
--anchor-len-deletion      290
--anchor-len-CGR           260
--anchor-len-duplication   210
--anchor-len-complexInsertion 230
--anchor-len-complexDeletion  290
--anchor-len-SCGR          340

--anchor-minDepth          0.33
--anchor-maxDepth          3.00
```

A mitochondrial chromosome or contig can optionally be excluded:

```bash
--mitochondria-chr <chromosome_name>
```

### Example with additional parameters

```bash
perl 01.extractSignature.pl \
    --bam sample.coordinate.sorted.bam \
    --bamSortRead sample.readname.sorted.bam \
    --output sample.signature \
    --min-mapq 20 \
    --min-sv-len 50 \
    --max-sv-len 10000
```

### SV signature classes

The current implementation includes the following classes:

```text
Insertion
Deletion
ComplexInsertion
ComplexDeletion
Duplication
CGR
SCGR
```

The output retains information describing:

```text
Chromosome
reference anchor coordinates
reference breakpoints
read name
read coordinates
read breakpoints
SV type
strand
anchor integrity
anchor scores
reference path
read path
original sequence
anchor/alternative sequence
```

These read-derived alternative sequences form the basis for subsequent pseudo-genome construction.

---

## Step 2. Cluster SV signatures

The extracted SV signatures are clustered using:

```text
02.cluster.pl
```

### Usage

```bash
perl 02.cluster.pl \
    --input sample.signature.gz \
    --output sample.cluster
```

The clustering step combines redundant SV signatures originating from different HiFi reads and identifies representative events.

### Default clustering distances

```text
Insertion      130 bp
Deletion        70 bp
CGR            360 bp
Duplication    220 bp
SCGR           410 bp
```

### Default minimum read-support thresholds

```text
Insertion        2
Deletion         2
CGR             13
Duplication     18
SCGR             18
```

These parameters can be adjusted using:

```text
--max-distance-insertion
--max-distance-deletion
--max-distance-CGR
--max-distance-duplication
--max-distance-SCGR

--min-reads-support-insertion
--min-reads-support-deletion
--min-reads-support-CGR
--min-reads-support-duplication
--min-reads-support-SCGR
```

Additional parameters include:

```text
--weight-complexINSDEL
--weight-unCover
--block-gap
--node-max-distance
--intermediate-output
```

### Example

```bash
perl 02.cluster.pl \
    --input sample.signature.gz \
    --output sample.cluster \
    --min-reads-support-insertion 2 \
    --min-reads-support-deletion 2
```

Principal output:

```text
sample.cluster.gz
```

The clustered output retains representative alternative sequences for subsequent pseudo-genome construction.

---

## Step 3. Generate pseudo-genomes

Pseudo-genomes are generated from the clustered SV events using:

```text
03.generatePseudogenome.pl
```

Index the reference genome first:

```bash
samtools faidx reference.fa
```

### Usage

```bash
perl 03.generatePseudogenome.pl \
    --sv_file sample.cluster.gz \
    --ref_genome reference.fa \
    --sv_distance 250 \
    --output_prefix sample
```

The default minimum distance between SV events assigned to the same pseudo-genome is:

```text
250 bp
```

PooledPP attempts to place compatible SV events within the same pseudo-genome.

If a candidate event cannot be placed into an existing pseudo-genome because it is too close to a previously assigned event, it is assigned to another pseudo-genome path.

The resulting pseudo-genomes are:

```text
sample.Pseudogenome1.fa
sample.Pseudogenome2.fa
sample.Pseudogenome3.fa
...
```

The corresponding SV information files are:

```text
sample.SV_info1.txt.gz
sample.SV_info2.txt.gz
sample.SV_info3.txt.gz
...
```

For each retained event, the representative read-derived sequence replaces the corresponding reference genomic interval.

The generated pseudo-genomes therefore encode compatible collections of alternative sequence paths recovered directly from pooled HiFi reads.

---

## Step 3.1. Construct the graph pangenome

`03.generatePseudogenome.pl` automatically generates a Minigraph command script:

```text
sample.minigraph.sh
```

The generated command follows the general form:

```bash
minigraph -cxggs -t20 \
    reference.fa \
    sample.PseudogenomeN.fa \
    ... \
    sample.Pseudogenome2.fa \
    sample.Pseudogenome1.fa \
    > sample.gfa
```

Run:

```bash
bash sample.minigraph.sh
```

Principal output:

```text
sample.gfa
```

The graph contains the reference genome backbone together with alternative sequence paths represented by the PooledPP pseudo-genomes.

---

## Step 4. Generate VCF from the graph

### ⚠ vg version limitation

This step has a strict vg compatibility requirement.

```text
Recommended:       vg v1.40.0
Compatible range:  vg <= v1.40.0
Not supported:     vg > v1.40.0
```

Before running this step, check the installed vg version:

```bash
vg version
```

The recommended version is:

```text
vg v1.40.0
```

### Run graph-to-VCF conversion

```bash
perl 04.generateVCF.pl \
    --input sample.gfa \
    --output sample.graph \
    --vg /path/to/vg
```

Equivalent short options are:

```bash
perl 04.generateVCF.pl \
    -i sample.gfa \
    -o sample.graph \
    -v /path/to/vg
```

The script first converts the GFA graph into VG format:

```bash
vg view -F -g sample.gfa > sample.graph.vg
```

It subsequently deconstructs the graph into VCF:

```bash
vg deconstruct sample.graph.vg -a > sample.graph.vcf
```

The generated VCF is further processed to:

1. add the `GT` FORMAT field;
2. assign the graph-derived alternative-path genotype;
3. compress and index the VCF;
4. normalize multiallelic variants using `bcftools norm`;
5. split multiallelic records into normalized bi-allelic representations; and
6. compress and index the final VCF.

Principal outputs include:

```text
sample.graph.addGT.vcf.gz
sample.graph.addGT.vcf.gz.tbi

sample.graph.split.vcf.gz
sample.graph.split.vcf.gz.tbi
```

> **Do not use vg versions newer than v1.40.0 for the current `04.generateVCF.pl` workflow.**

---

# Complete example

A complete PooledPP workflow is shown below.

## 1. Prepare BAM files

```bash
samtools sort \
    -o sample.coordinate.sorted.bam \
    sample.bam

samtools sort -n \
    -o sample.readname.sorted.bam \
    sample.bam

samtools index sample.coordinate.sorted.bam
```

## 2. Extract SV signatures

```bash
perl 01.extractSignature.pl \
    --bam sample.coordinate.sorted.bam \
    --bamSortRead sample.readname.sorted.bam \
    --output sample.signature
```

## 3. Cluster SV signatures

```bash
perl 02.cluster.pl \
    --input sample.signature.gz \
    --output sample.cluster
```

## 4. Generate pseudo-genomes

```bash
samtools faidx reference.fa

perl 03.generatePseudogenome.pl \
    --sv_file sample.cluster.gz \
    --ref_genome reference.fa \
    --sv_distance 250 \
    --output_prefix sample
```

## 5. Construct the graph pangenome

```bash
bash sample.minigraph.sh
```

## 6. Verify vg version

```bash
vg version
```

Recommended:

```text
vg v1.40.0
```

## 7. Generate graph-derived variants

```bash
perl 04.generateVCF.pl \
    --input sample.gfa \
    --output sample.graph \
    --vg /path/to/vg
```

---

## Output hierarchy

```text
Pooled HiFi alignment
        |
        v
sample.signature.gz
        |
        v
sample.cluster.gz
        |
        +-------------------------+
        |                         |
        v                         v
Pseudogenome1.fa          Pseudogenome2.fa ...
        |                         |
        +------------+------------+
                     |
                     v
                 Minigraph
                     |
                     v
                 sample.gfa
                     |
                     v
              vg <= v1.40.0
                     |
                     v
          sample.graph.vg
                     |
                     v
              vg deconstruct
                     |
                     v
         sample.graph.split.vcf.gz
```

---

## Methodological concept

The central idea of PooledPP is to convert **read-supported structural variation from pooled HiFi sequencing data into explicit sequence paths for graph pangenome construction**.

Instead of assigning pooled reads to specific individuals or attempting individual-level genome assembly, PooledPP:

1. extracts candidate structural variant paths from HiFi alignments;
2. evaluates local alignment and anchor information;
3. clusters SV signatures supported by multiple reads;
4. retains representative read-derived alternative sequences;
5. separates incompatible or closely positioned SV events among multiple pseudo-genomes;
6. inserts these alternative sequences into the corresponding reference genomic background;
7. constructs a graph from the reference genome and pseudo-genomes; and
8. deconstructs the resulting graph to obtain graph-derived variant representations.

The pseudo-genomes therefore serve as **computational carriers of alternative sequence paths** rather than reconstructed biological individuals or fully phased haplotypes.

---

## Important notes

### 1. Pseudo-genomes are not individual genomes

The pseudo-genomes generated by PooledPP do not correspond to specific sequenced individuals.

They represent computational combinations of compatible alternative sequence paths extracted from pooled HiFi reads.

### 2. PooledPP is not a haplotype assembler

PooledPP does not attempt to infer complete diploid genomes or chromosome-scale haplotypes from individual animals within the sequencing pool.

### 3. SV support thresholds depend on sequencing characteristics

The default clustering and read-support parameters were designed for the current PooledPP workflow. Users analyzing datasets with substantially different sequencing depth or pool size should evaluate whether parameter adjustment is appropriate.

### 4. vg compatibility is critical

`04.generateVCF.pl` was developed and tested using:

```text
vg v1.40.0
```

The current workflow should therefore be run with:

```text
vg <= v1.40.0
```

Versions newer than v1.40.0 are not supported by the current implementation because the graph-to-VCF conversion procedure required by PooledPP is no longer available in the same form in later vg releases.

For maximum reproducibility, **vg v1.40.0 is strongly recommended**.

### 5. Reproducibility

For publication-quality analyses, software versions and major parameters used for PooledPP should be recorded and reported together with the results.

---

## Citation

If you use PooledPP in your research, please cite the publication describing the method.

```text
Citation information will be added upon publication.
```

---

## License

Please refer to the `LICENSE` file for licensing information.

---

## Contributors

**Liqin Chen**  
**Qiuming Chen**

---

## Contact

For questions, bug reports, or feature requests, please open an issue in the GitHub repository.
