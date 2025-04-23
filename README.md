# Spermosphere_Across_Different_Soils
Analyses of the spermosphere microbiome from different soils across america.
Goal of this is to study if core microbes are present even across different geographic regions.

# Soil Collection
Soil was collected from: 
  -Southern Alabama (Wiregrass), 
  -Central Alabama (E.V. SMITH), 
  -Northern AL (TVREC), 
  -Michigan, 
  -Pennsylvania, 
  -Arkansas,
  -North Dakota,
  -Potting soil that was autoclaved thrice.
Information of the crop fields were recorded in the metadata file (Spermosphere_16S_Metadata.csv)

# Sampling
Soil was allowed to equilibrate with 1 mL of sterile water in the growth chamber for ~ 2 hours.
Seeds were weighed prior to planting.
Seeds were planting in soil (in a 12 well plate), and some wells did not have seeds planted, done in replicates of 6.
Seeds were allowed to germinate in the soil for 0, 8, 12, and 16 hours.
Once the time limit was reached, the seed and the spermosphere soil was taken out using a 11 mm cork borer and put into a sterile coin envelope.
Seeds were taken and weighed for total water imbibed.
Soil was colleccted into a 2 mL sterile tube and then frozen at -80 degrees C.
The soil was then extracted using the OmegaBiotek soil DNA extraction kit (Ivrine, CA).

# Sequencing
Both the V3 and V4 region were sequenced for the Bacterial 16S Sequences via Nextseq from SeqCenter LLC (Pittsburgh PA)

## Link to Analyses: These are links to .md files, other files are also available.
- [Bacteria Quality Trimming, Taxonomic Classification and OTU clustering]
- [Bacterial Sample Preproccessing](Bacteria/2_Bacteria_Preprocessing/LL_Prokaryote_Publication_Preprocessing.md)
- [Alpha and Beta Diversity Analysis](Bacteria/3_Alpha_and_Beta_Diversity/LL_Prokaryote_Diversity_Analysis.md)

# Bacterial Sequence Preprocessing (HPC)
First the in Output from the Bacteria quality was assessed via FastQC, 
reverse reads were merged via USEARCH (v12.0.beta)?, 
primers were stripped via CUTADAPTv210, 
reads were  then truncated at equal length and singletons and chimeras were removed and then clustered using USEARCH (v12.0.beta)?,
Prokaryote OTU’s were aligned with  MAFFT (v7.505) using the FFT-NS-2 algorithm and Fastree (v2.1.10 SSE3) was used to generate a phylogenetic tree,
Taxonomy was assigned to the pOTUs using the SINTAX algorithm with the SILVA 138.1 database.

Output Files: OTU table, FASTA file, Taxa Table, Phylogenetic Tree

# Bacterial Sample Preprocessing (R)
- [Bacterial Sample Preproccessing](Bacteria/2_Bacteria_Preprocessing/LL_Prokaryote_Publication_Preprocessing.md)


Output Files were read in as well as a metadata file into R.
NOTE: OTU Table and FASTA files are in a zip folder so if this project is downloaded this folder needs unzipped first and files need to be put into the directory with the other files.
Files were merged into a phlyoseq object

## Decontamination
Samples were decontaminated from the negative extraction controls and negative PCR controls via the package decontam_1.26.0.
Samples were filtered so that they only included bacterial taxa.

## Mock Community Analysis
To ensure that the illumina flow cell had minimal cross over occuring, the mock community analysis is used as a reassurance that data was not manipulated.
A mock community of known 8 members was sequenced along with the rest of the samples.
A theoretical community is also made with the ideal fractions of each bacteria.

## Rarefraction and Read-depth Analysis
This was done to assess sequencing depth and what was the median read depth across samples.

## Rarefication and CSS normalization
Samples were CSS Normalized and Rarefied to have the output files CSS normalized (Spermosphere_CSS.RDS) and Rarefied (Spermosphere_Rare.RDS)

# Alpha and Beta Diversity Analysis
- [Alpha and Beta Diversity Analysis](Bacteria/3_Alpha_and_Beta_Diversity/LL_Prokaryote_Diversity_Analysis.md)
## Alpha Diversity
Pielou's Evenness, Faith's Phylogenetic Diversity and Richness was calculated for each Location and Sample type.
For Richness statistcs a negative binomial model was used, for evenness a betareg linear model was used, for Faith's Phylogenetic Diversity, a negative binomial model was used.
Observed Richness, Evenness and Faith's Phylogenetic Diversity was also calculated
## Beta Diversity
Weighted and Unweighted Unifrac matrices were used for this analysis because it best answered the question of "What is the core"
Weighted Unifrac takes into account relative abundancies and phylogenetic information, whereas unweighted unifrac takes into account presence abscenses with the phylogenetic information.
A global PcOA plot, permanova was performed and since all interactions were significant, the data was subset based on time.
Betadispersion was calculated for each time point and type of Unifrac

#File Tree
```
.
├── Bacteria #Bacterial Analyses Processes
│   ├── 1_Bacteria_QualityTrimming_HPC #Processing reads by trimming, clustering and taxonomic classifications
│   │   ├── Bacteria_HPC
│   │   │   ├── HPC_Input
│   │   │   └── HPC_Output
│   │   │       ├── 16s_taxonomy_SINTAX.txt #OTU Taxonomic Classification
│   │   │       ├── FASTA _and_OTU_TABLE.zip #Together The Fasta (Sequences) and the OTU Table (OTU with reads per sample)
│   │   │       ├── otu_tree.tre #Phylogenetic tree for the sequences
│   │   │       └── Spermosphere_16S_Metadata.csv #Metadata
│   │   └── Bacteria_RDS
│   │       ├── physeq_bac_full.rds #RDS before cleaning 
│   │       ├── Spermosphere_bac_clean.rds #Cleaned and Decontaminated RDS
│   │       ├── Spermosphere_CSS.rds #CSS Normalized RDS
│   │       └── Spermosphere_Rare.rds #Rarefied RDS
│   ├── 2_Bacteria_Preprocessing
│   │   ├── Bacteria_RDS 
│   │   │   ├── physeq_bac_full.rds #RDS before cleaning
│   │   │   ├── Spermosphere_bac_clean.rds #Cleaned and Decontaminated RDS
│   │   │   ├── Spermosphere_CSS.rds #CSS Normalized RDS
│   │   │   └── Spermosphere_Rare.rds #Rarefieid RDS
│   │   ├── HPC_Output #put in this directory for good RMD file Pathing
│   │   │   ├── 16s_taxonomy_SINTAX.txt #See HPC output above
│   │   │   ├── FASTA _and_OTU_TABLE.zip #See HPC output above
│   │   │   ├── otu_tree.tre #See HPC output above
│   │   │   └── Spermosphere_16S_Metadata.csv #See HPC output above
│   │   ├── LL_Prokaryote_Publication_Preprocessing.html #HTML version of RMD
│   │   ├── LL_Prokaryote_Publication_Preprocessing.md #Github friendly version of RMD
│   │   ├── LL_Prokaryote_Publication_Preprocessing.Rmd #Code for the Preprocessing
│   │   ├── LL_Prokaryote_Publication_Preprocessing_files #Output from the md creation
│   │   │   └── figure-markdown_strict
│   │   │       ├── Mock Community Analysis-1.png
│   │   │       ├── Rarefaction Analysis-1.png
│   │   │       └── Rarefaction Analysis-2.png
│   │   └── Spermosphere_Rare.rds
│   └── 3_Alpha_and_Beta_Diversity
│       ├── Bacteria_RDS #put in this directory for good RMD file Pathing
│       │   ├── physeq_bac_full.rds #SEE ABOVE
│       │   ├── Spermosphere_bac_clean.rds #SEE ABOVE
│       │   ├── Spermosphere_CSS.rds #SEE ABOVE
│       │   └── Spermosphere_Rare.rds #SEE ABOVE
│       ├── LL_Prokaryote_Diversity_Analysis.html
│       ├── LL_Prokaryote_Diversity_Analysis.md
│       ├── LL_Prokaryote_Diversity_Analysis.Rmd #Code for the Alpha and Beta Diversity
│       └── LL_Prokaryote_Diversity_Analysis_files #Output from the md creation
│           ├── figure-html
│           │   ├── 0hours Weighted Unifrac-1.png
│           │   ├── 12hours Weighted Unifrac-1.png
│           │   ├── 16hours Weighted Unifrac-1.png
│           │   ├── 8hours Weighted Unifrac-1.png
│           │   ├── Global Weighted Unifrac-1.png
│           │   ├── Plotting Faiths Phylogenetic Diversity-1.png
│           │   ├── Plotting Faiths Phylogenetic Diversity-2.png
│           │   ├── Plotting Pielous Eveness-1.png
│           │   ├── Plotting Pielous Eveness-2.png
│           │   ├── Plotting Richness-1.png
│           │   └── Plotting Richness-2.png
│           └── figure-markdown_strict
│               ├── 0hours Unweighted Unifrac-1.png
│               ├── 0hours Weighted Unifrac-1.png
│               ├── 12hours Unweighted Unifrac-1.png
│               ├── 12hours Weighted Unifrac-1.png
│               ├── 16hours Unweighted Unifrac-1.png
│               ├── 16hours Weighted Unifrac-1.png
│               ├── 8hours Unweighted Unifrac-1.png
│               ├── 8hours Weighted Unifrac-1.png
│               ├── Global Unweighted Unifrac-1.png
│               ├── Global Weighted Unifrac-1.png
│               ├── Plotting Faiths Phylogenetic Diversity-1.png
│               ├── Plotting Faiths Phylogenetic Diversity-2.png
│               ├── Plotting Pielous Eveness-1.png
│               ├── Plotting Pielous Eveness-2.png
│               ├── Plotting Richness-1.png
│               └── Plotting Richness-2.png
├── LICENSE
├── README.md
├── Rmd #Future Analyses RMD
│   ├── LL_Assembly_Analysis.Rmd #Looking at assembly and how the core microbiome changes over time
│   ├── LL_Differntial_Abundance.Rmd #Looking at differentially abundant taxa and how it can relate to the core microbiome
│   ├── LL_Networking_Analysis.Rmd #Looking at hub taxa that may be keystone species within the microbial communities
│   └── Maps.Rmd  #Sampling Maps
├── Spermosphere_Across_Different_Soils.Rproj
└── Spermosphere_Metadata_030725.csv
```
#Packages Used
attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] multcompView_0.1-10  emmeans_1.11.0       betareg_3.2-2        MASS_7.3-65          picante_1.8.2       
 [6] nlme_3.1-164         ape_5.8-1            decontam_1.26.0      ggrepel_0.9.6        microbiome_1.28.0   
[11] Biostrings_2.74.1    GenomeInfoDb_1.42.3  XVector_0.46.0       IRanges_2.40.1       S4Vectors_0.44.0    
[16] ggpubr_0.6.0         metagenomeSeq_1.48.1 RColorBrewer_1.1-3   glmnet_4.1-8         Matrix_1.7-0        
[21] limma_3.62.2         Biobase_2.66.0       BiocGenerics_0.52.0  lubridate_1.9.4      forcats_1.0.0       
[26] stringr_1.5.1        dplyr_1.1.4          purrr_1.0.4          readr_2.1.5          tidyr_1.3.1         
[31] tibble_3.2.1         ggplot2_3.5.2        tidyverse_2.0.0      vegan_2.6-10         lattice_0.22-6      
[36] permute_0.9-7        phyloseq_1.50.0      BiocManager_1.30.25 

loaded via a namespace (and not attached):
 [1] Wrench_1.24.0           bitops_1.0-9            sandwich_3.1-1          rlang_1.1.6            
 [5] magrittr_2.0.3          multcomp_1.4-28         ade4_1.7-23             matrixStats_1.5.0      
 [9] compiler_4.4.1          flexmix_2.3-20          mgcv_1.9-1              vctrs_0.6.5            
[13] reshape2_1.4.4          pkgconfig_2.0.3         shape_1.4.6.1           crayon_1.5.3           
[17] fastmap_1.2.0           backports_1.5.0         caTools_1.18.3          rmarkdown_2.29         
[21] tzdb_0.5.0              UCSC.utils_1.2.0        modeltools_0.2-23       xfun_0.52              
[25] zlibbioc_1.52.0         jsonlite_2.0.0          biomformat_1.34.0       rhdf5filters_1.18.1    
[29] Rhdf5lib_1.28.0         broom_1.0.8             parallel_4.4.1          cluster_2.1.6          
[33] R6_2.6.1                stringi_1.8.7           car_3.1-3               estimability_1.5.1     
[37] lmtest_0.9-40           Rcpp_1.0.14             iterators_1.0.14        knitr_1.50             
[41] zoo_1.8-14              nnet_7.3-19             splines_4.4.1           igraph_2.1.4           
[45] timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8            
[49] yaml_2.3.10             gplots_3.2.0            codetools_0.2-20        plyr_1.8.9             
[53] withr_3.0.2             Rtsne_0.17              evaluate_1.0.3          survival_3.6-4         
[57] pillar_1.10.2           carData_3.0-5           KernSmooth_2.23-24      foreach_1.5.2          
[61] generics_0.1.3          hms_1.1.3               munsell_0.5.1           scales_1.3.0           
[65] gtools_3.9.5            glue_1.8.0              tools_4.4.1             data.table_1.17.0      
[69] locfit_1.5-9.12         ggsignif_0.6.4          mvtnorm_1.3-3           rhdf5_2.50.2           
[73] grid_4.4.1              colorspace_2.1-1        GenomeInfoDbData_1.2.13 Formula_1.2-5          
[77] cli_3.6.4               gtable_0.3.6            rstatix_0.7.2           digest_0.6.37          
[81] TH.data_1.1-3           htmltools_0.5.8.1       multtest_2.62.0         lifecycle_1.0.4        
[85] httr_1.4.7              statmod_1.5.0   

