# Phylloglossum RNA Editing
Open access data and methodology for detecting RNA editing events in the chloroplast and mitochondria of the Western Australian lycophyte Phylloglossum drummondii.

This repository contains the code used to create figures presented in The Plant Journal article "Insights into U-to-C editing from the lycophyte Phylloglossum drummondii", as well as code for identifying RNA editing sites within chloroplast and mitochondrion transcripts using the Julia program 'pyrimid', and associated RNA editing Jupyter notebooks which operate with a Julia kernal.

## Data Availability

## Chloroplast genome assembly and annotation
DNA reads were trimmed using BBDuk from the [BBtools suite](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (Settings ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo) (Bushnell, 2016). Trimmed reads were used to assemble the chloroplast genome using NOVOPlasty with the Huperzia serrata rbcL gene as a seed sequence (settings type=chloro; genome range=120000-220000bp; kmer=51) (Dierckxsens et al., 2017), and using GetOrganelle (Jin et al., 2020). The resulting assemblies were identical. The chloroplast genome assembly was verified with Pilon (Walker et al., 2014). Using the RNA editing sites identified in the chloroplast transcriptome, we created an edited copy of the plastome where C was modified to T at sites exceeding 20% editing. This ensured that start and stop codons created by editing were accounted for. Annotation of the chloroplast genes was then done by Chloë (https://github.com/ian-small/chloe) and applied to the ‘unedited’ original version of the plastome. The chloroplast genome of 
Phylloglossum drummondii is available in GenBank as accession OR992133.

## Mitochondrial genome assembly and annotation
Mitochondrial assembly was achieved using a combined approach utilising NOVOPlasty (Dierckxsens et al., 2017), GetOrganelle (Jin et al., 2020), SPAdes (Bankevich et al., 2012) and Geneious Prime (v2023.0.4) (https://www.geneious.com). NOVOPlasty was able to obtain a non-contiguous mitochondrial assembly (settings type=mito_plant; genome range=300000-550000bp; kmer=23) using the complete mitochondrial genome of Phlegmariurus squarrosus as a seed input file. The NOVOPlasty assembly of the Phylloglossum chloroplast genome was supplied as a reference sequence. SPAdes  was used to assemble an alternative view of the mitochondrial genome as an assembly graph (settings --cov-cutoff 50.0 --assembler-only --careful) (Bankevich et al., 2012). The SPAdes output FastG graph was visualised in Bandage (Wick et al., 2015) and mitochondrial genes were identified by blastn (Altschul et al., 1990) using gene sequences from Phlegmariurus squarrosus as queries. DNA reads were mapped to each mitochondrial assembly option and the mapped reads were used as input for the Geneious de novo assembly algorithm. Ultimately, the SPAdes assembly version was used as the basis for the final assembly using connections suggested by Geneious and/or NOVOPlasty to manually edit the FastG file and eliminate dead-ends. Pilon (Walker et al., 2014) was used to verify the final assembly. The path chosen for the ‘master circle’ view of the mitochondrial genome is only one of many possible arrangements of the assembly. The Phylloglossum drummondii mitochondrial genome is available in GenBank as accession PP024676.

Annotations of the mitochondrial genome of Phlegmariurus squarrosus (GenBank accession JQ002659) were extracted using Geneious Prime (v2023.0.4) (https://www.geneious.com) and mapped to the Phylloglossum mitochondrial genome assembly. Intron-containing genes were problematic due to varying intron length between Phlegmariurus squarrosus and Phylloglossum. In these cases, individual exons were mapped and intron boundaries were verified using mapped RNA reads. tRNAscan-SE (Chan and Lowe, 2019) was used to check for missed tRNA gene annotations and tRNAs were checked against the PlantRNA2.0 database (Cognat et al., 2022). Genes were manually curated with reference to RNA editing events to ensure start codon creation, stop codon creation and premature stop codon removal events were taken into account. Annotations from Phylloglossum were mapped to the mitochondrial genomes of Huperzia crispata and Phlegmariurus squarrosus to compare annotation start and end points, and to compare gene counts.

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
