PooledPP enables accurate, fast and low-cost construction of pseudo-pangenome graphs with pooled PacBio HiFi data

![Figure 1](https://github.com/user-attachments/assets/1966db69-6de1-40cf-b354-29ff5b2dd42e)
Workflow of PooledPP

PooledPP is a robust pipeline for constructing a pseudo-pangenome from pooled HiFi sequencing data.  
The workflow begins by aligning HiFi reads to a reference genome using pbmm2. Structural variation (SV) signatures are extracted (01.extractSignature.pl) and clustered (02.cluster.pl). Clustered signatures are then used to generate pseudo-genomes (03.generatePseudogenome.pl), which are further assembled into a pseudo-pangenome graph using minigraph (panel a).  
For uniquely aligned reads, PooledPP detects two canonical SV types:  
•	Deletions: regions in the reference that are skipped by the read alignment (panel b);  
•	Insertions: sequences present in the read but absent from the reference (panel c).  
For split-read alignments (reads with two or more alignment segments), SVs are further categorized based on their mapping context:  
•	Complex deletions and complex insertions follow the same biological definitions as simple deletions and insertions, but are located within split alignments (panels d and e).  
When two adjacent segments of a read are separated by both a large insertion and a large deletion, the event is classified as a complex genomic rearrangement (CGR) (panel f), representing more extensive and compound structural changes.  
Duplications are identified as split-read events where adjacent alignment segments on the read overlap in the reference coordinate, indicating that a segment of the reference is repeated. These duplications may represent tandem or interspersed duplications.  
If a read contains three or more alignment segments, and any pair of non-adjacent segments (e.g., segments 1 and 3 or 2 and 5) show consistent strand orientation and continuity across both reference and read coordinates, the event is classified as a super complex genomic rearrangement (SCGR) (panel g). These signatures capture higher-order SVs including inversions and translocations.  
To refine SV breakpoints, PooledPP employs an anchor-based strategy (panel h): For each breakpoint, local anchor sequences with 100% identity are searched bidirectionally. Anchors are scored based on their perfect match length and local read depth, normalized against genome-wide depth, to prioritize high-confidence breakpoints.  
For duplications and SCGRs, PooledPP represents alignments using a graph-based model, where nodes correspond to alignment segments and edges denote continuity in both read and reference coordinates (panel i). This enables robust resolution of complex rearrangements in repetitive or mosaic genomic regions.  
