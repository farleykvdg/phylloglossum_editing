# _Phylloglossum_ RNA Editing
Open access data and methodology for detecting RNA editing events in the chloroplast and mitochondria of the Western Australian lycophyte Phylloglossum drummondii.

This repository contains the code used to create figures presented in The Plant Journal article "Insights into U-to-C editing from the lycophyte Phylloglossum drummondii", as well as code for identifying RNA editing sites within chloroplast and mitochondrion transcripts using the Julia program 'pyrimid', and associated RNA editing Jupyter notebooks which operate with a Julia kernal.

## Data Availability

DNAseq and RNAseq datasets are available on the Sequence Read Archive as accession [PRJNA818771](https://www.ncbi.nlm.nih.gov/sra?term=SRP365360).

The chloroplast genome of _Phylloglossum drummondii_ is available in GenBank as accession [OR992133](https://www.ncbi.nlm.nih.gov/nuccore/OR992133).

The mitochondrion genome of _Phylloglossum drummondii_ is available in GenBank as accession [PP024676](https://www.ncbi.nlm.nih.gov/nuccore/PP024676).

## Chloroplast genome assembly and annotation
Trim DNA reads using BBDuk from the [BBtools suite](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/).

Settings: _ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo_ 

Use trimmed reads to assemble the _Phylloglossum drummondii_ chloroplast genome using [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty). Use a related species chloroplast sequence as a seed, for example, [_Huperzia serrata rbcL_](https://www.ncbi.nlm.nih.gov/nuccore/DQ464224.1). 

Config file settings: _type=chloro; genome range=120000-220000bp; kmer=51_

Optionally, you can verify and attempt to improve the assembly using [Pilon](https://github.com/broadinstitute/pilon). 

Prelimintary annotation of the chloroplast genome can be achieved using [Chloë](https://chloe.plastid.org/annotate.html), however, Chloë is optimised for angiosperms and will not currently produce accurate annotations. Annotations for the _Phylloglossum drummondii_ plastome were generated with a beta version of Chloë which included Lycophyte and Monilophyte reference sequences.

## Mitochondrial genome assembly and annotation
Mitochondrial assembly can be achieved by combining NOVOPlasty, [SPAdes](https://github.com/ablab/spades) and [Geneious Prime](https://www.geneious.com/). 

NOVOPlasty is capable of assembling non-contiguous mitochondrial assemblies using Settings _type=mito_plant; genome range=300000-550000bp; kmer=23_. We supplied the complete mitochondrial genome of [_Phlegmariurus squarrosus_](https://www.ncbi.nlm.nih.gov/nucleotide/NC_017755.1) as a seed input file, and our own assembly of the _Phylloglossum drummondii_ plastome as a secondary reference in the NOVOPlasty .config file. 

SPAdes was used to assemble the mitochondrial genome using Settings _--cov-cutoff 50.0 --assembler-only --careful_.

The SPAdes output FastG graph was visualised in [Bandage](https://github.com/rrwick/Bandage) and mitochondrial genes were identified by blastn using sequences from _Phlegmariurus squarrosus_ as queries. 

DNA reads were mapped to each mitochondrial assembly option and the mapped reads were used as input for the Geneious _de novo_ assembly algorithm. Ultimately, the SPAdes assembly version was used as the basis for the final assembly using connections suggested by Geneious and/or NOVOPlasty to manually edit the FastG file and eliminate dead-ends. Pilon was used to verify the final assembly. The path chosen for the ‘master circle’ view of the mitochondrial genome is only one of many possible arrangements of the assembly.

Annotations from [_Phlegmariurus squarrosus_](https://www.ncbi.nlm.nih.gov/nucleotide/NC_017755.1) were extracted using Geneious Prime and mapped to the Phylloglossum mitochondrial genome assembly. You can also use Geneious' "annotate from" function to quickly map gene annotations from related species onto _de novo_ assemblies. In cases where introns contained large insertions or deletions, mapping individual exons accurately identified intron boundaries, which were later verified with RNAseq data. 

[tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) was used to check for missed tRNA gene annotations, and identified tRNAs were then checked against the [PlantRNA2.0 database](https://seve.ibmp.unistra.fr/plantrna/). 

All mitochondrion genes were manually curated using RNA editing events as a reference to ensure start codon creation, stop codon creation and premature stop codon removal events were accounted for. 

We subsequently mapped annotations from _Phylloglossum_ to the mitochondrion genomes of _Huperzia crispata_ and _Phlegmariurus squarrosus_ to compare annotation start points, end points and gene counts. We have corrected several annotations in the _Huperzia_ and _Phlegmariurus_ mitochondrion genomes.

## Detection of RNA editing events
RNA editing events were determined by mapping trimmed DNAseq and RNA-seq data to assembled chloroplast and mitochondrial genomes by merging forward and reverse read files with BBmerge (settings qtrim2=t, trimq=10,15,20, minq=12) and mapping merged and unmerged reads to the genome assemblies using BBWrap (settings mappedonly=t; ambiguous=random) (Bushnell, 2016). Nucleotide count files for each position in the organelle genomes were generated for each RNA-seq dataset and the DNAseq dataset using Pyrimid (settings –m 0, -u) (https://github.com/ian-small/pyrimid). DNA and RNA base identities at each position of the chloroplast and mitochondrial genomes were filtered for pyrimidine mismatches. Ribosomal RNA genes were masked as the RNA samples were depleted of plant rRNAs but may contain contaminating rRNA from microorganisms. For mitochondrial reads, binomial tests were performed on the mapped and unmapped RNA reads to calculate whether editing significantly exceeded an arbitrary threshold of 5% at each site. P-values were corrected for multiple testing using the Benjamini-Hochberg procedure. For chloroplast reads, Fisher’s exact test was also used at each site to determine whether the proportion of ‘edited’ reads was significantly greater in the RNA-seq data than in the DNA-seq data. For the mitochondrial reads, this test was not used as the lower coverage of DNA reads mapping to the mitochondrial genome made the test insensitive. A full list of RNA editing events identified is available in Supplementary Table S2. RNA editing events are annotated as misc_features of both the chloroplast and mitochondrial genome GenBank entries. 
De novo transcriptome assembly
Using RNA reads trimmed with BBDuk (settings ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo) (Bushnell, 2016), de novo transcriptome assemblies for Phylloglossum, Huperzia serrata and Phlegmariurus squarrosus were built using rnaSPAdes run with default settings (Bushmanova et al., 2019). RNA-seq data for Huperzia serrata was downloaded from Genome Sequence Archive accession number PRJCA000351 (Yang et al., 2017). RNA-seq data for Phlegmariurus squarrosus was downloaded from medplantRNA-seq.org (https://medplantrnaseq.org/assemblies/Huperzia_squarrosa.tar.gz).

## Identification of PPR proteins
Open reading frames in the Phylloglossum, Phlegmariurus squarrosus and Huperzia serrata transcriptomes were translated in forward and reverse orientations and in all six reading frames using orfinder.jl (github.com/ian-small/orfinder). “hmmsearch” from HMMER v3.2.1 (hmmer.org) was used to search for PPR motifs in the ORFs of the Phylloglossum transcriptome using the DYW and DYW:KP motif Hidden Markov Model (HMM) profiles from (Gutmann et al., 2020). PPRfinder (Gutmann et al., 2020) was used to predict the motif structure of each protein from the output of hmmsearch. 

## Phylogenetic analysis of DYW:KP proteins
Four DYW:KP sequences from Huperzia serrata, four sequences from Phylloglossum and three sequences from Phlegmariurus squarrosus and were extracted from assembled transcriptomes. The 11 DYW:KP domain peptide sequences were aligned with MUSCLE v5.1 (Edgar, 2004) and then converted to codon alignments. A maximum likelihood phylogenetic tree was constructed using IQTREE 2 v2.1.4 (Minh et al., 2020). The MG+F3X4+R2 model was identified as the best fit according to the inbuilt ModelFinder (Kalyaanamoorthy et al., 2017). Node confidence was estimated from 1000 bootstrap replicates with ultrafast bootstrap approximation UFBoot2 (Hoang et al., 2018).

## U-to-C PPR editing factor target prediction
The longest KP1, KP2, KP3 and KP4 sequences from Phylloglossum, Huperzia serrata and Phlegmariurus squarrosus were identified and supplied alongside the nucleotide sequences 30 nt upstream of mitochondrial U-to-C editing sites as input for the program PPRmatcher (Royan et al., 2021) which assesses the strength of a PPR motif match to a nucleotide using a motif scoring table derived from verified PPR binding sites. Proteins from each species were matched against the editing site from their respective mitochondrial genes, which were identical except for the nad5 upstream region. In some cases, there were not enough PPR motifs to make strong predictions, and the association of truncated DYW:KP proteins was made through the protein alignments, and codon model phylogeny, such as in the case of Phylloglossum KP1.

## Acknowledgements
This work was supported by Australian Research Council grant DP200102981. The authors have no conflicts of interest to declare. We are grateful to Volker Knoop and Simon Zumkeller for discussions and advice concerning lycophyte organelle gene structure and splicing.
