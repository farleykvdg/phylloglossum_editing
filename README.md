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

Annotations from _Phlegmariurus squarrosus_ were extracted using Geneious Prime and mapped to the Phylloglossum mitochondrial genome assembly. You can also use Geneious' "annotate from" function to quickly map gene annotations from related species onto _de novo_ assemblies. In cases where introns contained large insertions or deletions, mapping individual exons accurately identified intron boundaries, which were later verified with RNAseq data. 

[tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) was used to check for missed tRNA gene annotations, and identified tRNAs were then checked against the [PlantRNA2.0 database](https://seve.ibmp.unistra.fr/plantrna/). 

All mitochondrion genes were manually curated using RNA editing events as a reference to ensure start codon creation, stop codon creation and premature stop codon removal events were accounted for. 

We subsequently mapped annotations from _Phylloglossum_ to the mitochondrion genomes of _Huperzia crispata_ and _Phlegmariurus squarrosus_ to compare annotation start points, end points and gene counts. We have corrected several annotations in the _Huperzia_ and _Phlegmariurus_ mitochondrion genomes.

## Detection of RNA editing events
Merge trimmed RNAseq datasets using BBmerage with Settings _qtrim2=t, trimq=10,15,20, minq=12_ 

Map the merged _and_ unmerged reads to the organelle genome assemblies using BBWrap Settings _mappedonly=t ambiguous=random_. 

Nucleotide count files for each position in the organelle genomes were generated for each RNA-seq dataset and the DNAseq dataset using the version of Pyrimid available through this repository using settings _–m 0, -u_. 

Nucleotide count files can be converted to RNA editing tables using the Jupyter notebooks present in this repository. In short, the notebooks require .fasta and .gff files for the organelle assemblies, and the output from pyrimid (tab separated nucleotide count file). The workflow presented in the Jupyter notebooks require a [Julia](https://julialang.org/) kernel. To install this, follow instructions for [installing IJulia](https://www.geeksforgeeks.org/add-julia-kernel-to-jupyter/). The notebooks will mask rRNA and tRNAs, and perform binomial and Fisher's exact test on pyrimidine mismatches within the nucleotide count file. Sites passing these tests will be called as RNA editing events, and the statistics will be summarised in the resulting  .tsv file.

## De novo transcriptome assembly

Use trimmed RNA reads to build a _de novo_ transcriptome assembly for _Phylloglossum drummondii_ using rnaSPAdes run with default settings. 

We assembled transcriptomes for _Phylloglossum drummondii_, and also _Phlegmariurus squarrosus_ and _Huperzia serrata_. RNA-seq data for _Huperzia serrata_ is available under accession number [PRJCA000351](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA00035) and _Phlegmariurus squarrosus_ is available [here](https://medplantrnaseq.org/assemblies/Huperzia_squarrosa.tar.gz).

## Identification of PPR proteins
Open reading frames in the _Phylloglossum drummondii_, _Phlegmariurus squarrosus_ and _Huperzia serrata_ transcriptomes were translated in forward and reverse orientations and in all six reading frames using the version of orfinder.jl available in this repository. 

Use “hmmsearch” from [HMMER v3.2.1](hmmer.org) was to search for PPR motifs in the generated open reading frames .fasta file from the transcriptome assemblies using the DYW and DYW:KP motif Hidden Markov Model (HMM) profile "all_KP.hmm" present in this repository, and originally described in [(Gutmann et al., 2020)](https://doi.org/10.1016/j.molp.2019.11.002). 

Use PPRfinder_vApr19 (present in this repository) to identify PPR proteins from the open reading frames with the output of hmmsearch.

## Phylogenetic analysis of DYW:KP proteins
Four DYW:KP sequences from Huperzia serrata, four sequences from Phylloglossum and three sequences from Phlegmariurus squarrosus and were extracted from assembled transcriptomes. The 11 DYW:KP domain peptide sequences were aligned with MUSCLE v5.1 (Edgar, 2004) and then converted to codon alignments. A maximum likelihood phylogenetic tree was constructed using IQTREE 2 v2.1.4 (Minh et al., 2020). The MG+F3X4+R2 model was identified as the best fit according to the inbuilt ModelFinder (Kalyaanamoorthy et al., 2017). Node confidence was estimated from 1000 bootstrap replicates with ultrafast bootstrap approximation UFBoot2 (Hoang et al., 2018).

## U-to-C PPR editing factor target prediction
The longest KP1, KP2, KP3 and KP4 sequences from Phylloglossum, Huperzia serrata and Phlegmariurus squarrosus were identified and supplied alongside the nucleotide sequences 30 nt upstream of mitochondrial U-to-C editing sites as input for the program PPRmatcher (Royan et al., 2021) which assesses the strength of a PPR motif match to a nucleotide using a motif scoring table derived from verified PPR binding sites. Proteins from each species were matched against the editing site from their respective mitochondrial genes, which were identical except for the nad5 upstream region. In some cases, there were not enough PPR motifs to make strong predictions, and the association of truncated DYW:KP proteins was made through the protein alignments, and codon model phylogeny, such as in the case of Phylloglossum KP1.

## Acknowledgements
This work was supported by Australian Research Council grant DP200102981. The authors have no conflicts of interest to declare. We are grateful to Volker Knoop and Simon Zumkeller for discussions and advice concerning lycophyte organelle gene structure and splicing.
