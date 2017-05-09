---
layout: page
title: CBB752 Spring 2017
tagline: Final Project
---

2.1 Identifying off target CRISPR sites
------------------


Table of Contents
-----------------------




**Contributors**
 -Writing: Ramya
 -Coding: Jay
 -Pipeline: Daniel

### Introduction: CRISPR/Cas in the Context of Genome Editing
Since the discovery of its versatile genome editing function in 2013, the CRISPR/Cas system has been heavily studied for its use in targeted genome modification in eukaryotic cells. In many current applications, the Cas9 endonuclease from Streptococcus pyogenes is directed by a synthetic single stranded guide RNA (gRNA) containing homology to a genomic region. The gRNA displaces the duplex DNA and forms its own RNA:DNA duplex with one of the DNA strands. Cas9 then cleaves both strands of the DNA, inducing native DNA repair pathways. At this point, if a retinue of mutations at the site is the end goal, non-homologous end joining occurs with indels. However, if a specific point mutation or sequence insert is the goal, the addition of excess double stranded DNA of interest will result in a small population of recombinants with the change after Cas9 action.
< CRISPR image from addgene (show with HR)>

### Writing:

#### Off-Target Effects and “Determinants of CRISPR/Cas sgRNA Specificity” (1)
The sheer size of the human genome often means that the gRNA target sequence is at least partially homologous to several locations in the genome. This binding of gRNAs to unintended genomic sequences is referred to as off-target effects. Even with the careful design of gRNAs, off-target binding is high and gRNA dependent – between 10-1300 off target binding events were observed in one study with 12 different gRNAs (2). Currently, off-target effects are one of the primary set backs for clinical application of CRISPR/Cas9. However, there are several mechanisms that confer specificity to the CRISPR/Cas9 system that can be leveraged to reduce off-target effects.

*PAM site:*

As shown in orange in the figure above, the protospacer adjacent motif, or PAM site, is a 3 base pair sequence expected to bind the 3’ end of the designed gRNA. It has been shown to be absolutely required in the genome to target Cas9, as it is at the 5’ end of the PAM site that Cas9 introduces a double stranded break into the DNA. About 50% of the time in humans, the PAM sequence is the GGG motif, though NGG is also seen (2).  Depending on the origin species of the Cas9 employed, different PAM motifs can be probed. 

*gRNA: PAM Proximal versus PAM Distal*

The gRNA consists of two main elements: a scaffolding piece to bind Cas9 and a target piece complementary to a region in the genome. The designed gRNA is often about 20 nucleotides long. While intuition suggests that longer sequences may confer higher specificity, this was proven wrong (1). Free single stranded RNA is highly unstable, and it is thought that Cas9 protects only about 20 nucleotides from degradation. 

The PAM proximal region of gRNA is that which is closest to the PAM site while the PAM distal gRNA region is farther from the site. In a genome-wide ChIP-seq experiment using 12 gRNAs to characterize their off-target binding preferences, near-perfect homology of gRNA to DNA was observed at the PAM-proximal region (2). On the other hand, the PAM distal region had up to 10 mismatches with the genomic DNA (2). The PAM proximal region is referred to as the seed sequence to ensure gRNA specificity. The exact number of bases for this region is heavily debated but ranges between 1-12 nucleotides (1)(2)(3). 

The gRNA is the most important determinant of Cas9 specificity. The exact mechanism by which gRNA confers this specificity is not well understood. Some potential mechanisms include changing “the thermodynamic stability of the gRNA:DNA duplex,” the effective concentration of Cas9-gRNA (higher = less specific), the blocking of other Cas9 sites for trans-acting binding proteins, and the conformation of Cas9 to better access chromatin (1). 

*Chromatin Structure and Methylation:*

DNA exists as tightly compacted coils around histone proteins in human cells, and this results in bulky, complex structures that can sometimes impede protein-DNA interaction. Thus, gRNA binding is heavily biased towards binding more accessible regions of the genome. One study showed that over 50% of gRNA binding is in open chromatin regions (2). Less accessible chromatin regions also correspond to higher CpG methylation. The presence of these marks has been suggested to block Cas9 binding (1)(3).
*Cas9 related factors:*

Because of the fast action of Cas9 and rapid degradation, the addition of purified Cas9 is reported to help reduce off-target effects (3). Furthermore, Cas9 can be modified into a nickase, such that it only cuts one DNA strand and a pair of gRNAs and nickases can be used for greater stringency in mutation (4).


#### The Inputs Required for Off-Target-Effect analysis of gRNAs
For off-target analysis, most tools simply require the gRNA sequence and a genome to search. gRNA specific parameters can either be inputs or set by the program such as tolerated mismatches, length, etc. Some tools also allow the specification of PAM site based on the Cas9 ortholog used. More advanced tools allow for the inclusion of ENCODE data, specifying chromatin structure across the genome to account for this as well.

#### Comparison of De novo tool and Other Prediction tools such as: (Daniels)



### Coding:


#### Documentation: *De novo* Off-Target Mutation Prediction Tool for SubjectZ
The following is a basic overview of Jay's code.

#### Results:


### Pipeline:


#### Documentation:


#### Results:









#### Conclusions:








#### References:

 References can be included here or at the end of each relevant section.
 
 
