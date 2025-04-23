-   [Set global options](#set-global-options)
-   [Alpha Diversity Analysis](#alpha-diversity-analysis)
    -   [Autoclaved Soil](#autoclaved-soil)
    -   [EVSMITH Soil](#evsmith-soil)
    -   [Michigan Soil](#michigan-soil)
    -   [North Dakota Soil](#north-dakota-soil)
    -   [Wiregrass Soil](#wiregrass-soil)
    -   [Pennsylvania Soil](#pennsylvania-soil)
    -   [Arkansas Soil](#arkansas-soil)
    -   [TVREC Soil](#tvrec-soil)
    -   [Faith’s Phylogenetic Diversity
        Statistics](#faiths-phylogenetic-diversity-statistics)
-   [Beta Diversity](#beta-diversity)
    -   [Beta-Diversity using Weighted Unifrac Only
        Bacteria](#beta-diversity-using-weighted-unifrac-only-bacteria)
        -   [8 hours; Unifrac](#hours-unifrac)
        -   [12; Unifrac](#unifrac)
        -   [16; Unifrac](#unifrac-1)
    -   [Beta-Diversity using UnWeighted Unifrac Only
        Bacteria](#beta-diversity-using-unweighted-unifrac-only-bacteria)
        -   [8 hours; Unweighted Unifrac](#hours-unweighted-unifrac)
        -   [12; Unweighted Unifrac](#unweighted-unifrac)
        -   [16; Unweighted Unifrac](#unweighted-unifrac-1)
    -   [combining beta diversity
        figures](#combining-beta-diversity-figures)

\#Data Preprocessing# \## Libraries \##

    library(BiocManager)
    library(phyloseq)
    library(vegan)
    library(tidyverse)
    library(metagenomeSeq)
    library(ggpubr)
    library(Biostrings)
    library(microbiome)
    library(ggrepel)
    library(decontam)
    library(picante) #For Faiths Phylogentic Diversity
    library(MASS) #For Linear Modeling
    library(betareg) #For Linear Modeling
    library(emmeans)

## Set global options

    # Set global options #
    # no scientific notation
    options(scipen=10000) 

    # color blind pallets used throughout 
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
    tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

\##NOTE: Rmd files set the working directory to the place where the rmd
file is, not the actual working directory \##Read in and Separation##
\## Read in the RDS files \#

    ###Bacterial RDS
    BSpermosphere_RDS_CSS <- readRDS("Bacteria_RDS/Spermosphere_CSS.rds") ##Using Normalized Reads

    BSpermosphere_RDS_NoNorm <- readRDS("Bacteria_RDS/Spermosphere_bac_clean.rds") ##Using Non-normalized Reads

\##Seperation by Location \##

    bac.Auto <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "Autoclaved_Potting_Soil")
    bac.EV <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "EVSmith_EC")
    bac.MI <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "Michigan")
    bac.ND <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "N_Dakota")
    bac.Wire <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "Wiregrass_EC")
    bac.PA <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "Pennsylvania")
    bac.AR <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "Arkansas")
    bac.TVR <- subset_samples(BSpermosphere_RDS_NoNorm, Location_EC == "TVRec_EC")

# Alpha Diversity Analysis

## Autoclaved Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.Auto@sam_data$shannon <- estimate_richness(bac.Auto, measures=c("Shannon"))$Shannon
    bac.Auto@sam_data$invsimpson <- estimate_richness(bac.Auto, measures=c("InvSimpson"))$InvSimpson
    bac.Auto@sam_data$simpson <- estimate_richness(bac.Auto, measures=c("Simpson"))$Simpson
    bac.Auto@sam_data$richness <- estimate_richness(bac.Auto, measures=c("Observed"))$Observed
    bac.Auto@sam_data$even <- bac.Auto@sam_data$shannon/log(bac.Auto@sam_data$richness)
    Autosample.data.bac <- data.frame(bac.Auto@sam_data)
    Autosample.data.bac$Time <- factor(Autosample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    Auto_otu <- bac.Auto@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    autophylo.diversity <- pd(t(Auto_otu), bac.Auto@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    autophylo.diversity$Samples <- rownames(autophylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    autophylo.diversity.join <- left_join(autophylo.diversity, as.data.frame(bac.Auto@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    autophylo.diversity.join$Time.Point <- factor(autophylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## EVSMITH Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###
    bac.EV@sam_data$shannon <- estimate_richness(bac.EV, measures=c("Shannon"))$Shannon
    bac.EV@sam_data$invsimpson <- estimate_richness(bac.EV, measures=c("InvSimpson"))$InvSimpson
    bac.EV@sam_data$simpson <- estimate_richness(bac.EV, measures=c("Simpson"))$Simpson
    bac.EV@sam_data$richness <- estimate_richness(bac.EV, measures=c("Observed"))$Observed
    bac.EV@sam_data$even <- bac.EV@sam_data$shannon/log(bac.EV@sam_data$richness)
    EVsample.data.bac <- data.frame(bac.EV@sam_data)
    EVsample.data.bac$Time <- factor(EVsample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    EV_otu <- bac.EV@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    EVphylo.diversity <- pd(t(EV_otu), bac.EV@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    EVphylo.diversity$Samples <- rownames(EVphylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    EVphylo.diversity.join <- left_join(EVphylo.diversity, as.data.frame(bac.EV@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    EVphylo.diversity.join$Time.Point <- factor(EVphylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## Michigan Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.MI@sam_data$shannon <- estimate_richness(bac.MI, measures=c("Shannon"))$Shannon
    bac.MI@sam_data$invsimpson <- estimate_richness(bac.MI, measures=c("InvSimpson"))$InvSimpson
    bac.MI@sam_data$simpson <- estimate_richness(bac.MI, measures=c("Simpson"))$Simpson
    bac.MI@sam_data$richness <- estimate_richness(bac.MI, measures=c("Observed"))$Observed
    bac.MI@sam_data$even <- bac.MI@sam_data$shannon/log(bac.MI@sam_data$richness)
    MIsample.data.bac <- data.frame(bac.MI@sam_data)
    MIsample.data.bac$Time <- factor(MIsample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    MI_otu <- bac.MI@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    MIphylo.diversity <- pd(t(MI_otu), bac.MI@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    MIphylo.diversity$Samples <- rownames(MIphylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    MIphylo.diversity.join <- left_join(MIphylo.diversity, as.data.frame(bac.MI@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    MIphylo.diversity.join$Time.Point <- factor(MIphylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## North Dakota Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.ND@sam_data$shannon <- estimate_richness(bac.ND, measures=c("Shannon"))$Shannon
    bac.ND@sam_data$invsimpson <- estimate_richness(bac.ND, measures=c("InvSimpson"))$InvSimpson
    bac.ND@sam_data$simpson <- estimate_richness(bac.ND, measures=c("Simpson"))$Simpson
    bac.ND@sam_data$richness <- estimate_richness(bac.ND, measures=c("Observed"))$Observed
    bac.ND@sam_data$even <- bac.ND@sam_data$shannon/log(bac.ND@sam_data$richness)
    NDsample.data.bac <- data.frame(bac.ND@sam_data)
    NDsample.data.bac$Time <- factor(NDsample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    ND_otu <- bac.ND@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    NDphylo.diversity <- pd(t(ND_otu), bac.ND@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    NDphylo.diversity$Samples <- rownames(NDphylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    NDphylo.diversity.join <- left_join(NDphylo.diversity, as.data.frame(bac.ND@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    NDphylo.diversity.join$Time.Point <- factor(NDphylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## Wiregrass Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.Wire@sam_data$shannon <- estimate_richness(bac.Wire, measures=c("Shannon"))$Shannon
    bac.Wire@sam_data$invsimpson <- estimate_richness(bac.Wire, measures=c("InvSimpson"))$InvSimpson
    bac.Wire@sam_data$simpson <- estimate_richness(bac.Wire, measures=c("Simpson"))$Simpson
    bac.Wire@sam_data$richness <- estimate_richness(bac.Wire, measures=c("Observed"))$Observed
    bac.Wire@sam_data$even <- bac.Wire@sam_data$shannon/log(bac.Wire@sam_data$richness)
    Wiresample.data.bac <- data.frame(bac.Wire@sam_data)
    Wiresample.data.bac$Time <- factor(Wiresample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    Wire_otu <- bac.Wire@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    Wirephylo.diversity <- pd(t(Wire_otu), bac.Wire@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    Wirephylo.diversity$Samples <- rownames(Wirephylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    Wirephylo.diversity.join <- left_join(Wirephylo.diversity, as.data.frame(bac.Wire@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    Wirephylo.diversity.join$Time.Point <- factor(Wirephylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## Pennsylvania Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.PA@sam_data$shannon <- estimate_richness(bac.PA, measures=c("Shannon"))$Shannon
    bac.PA@sam_data$invsimpson <- estimate_richness(bac.PA, measures=c("InvSimpson"))$InvSimpson
    bac.PA@sam_data$simpson <- estimate_richness(bac.PA, measures=c("Simpson"))$Simpson
    bac.PA@sam_data$richness <- estimate_richness(bac.PA, measures=c("Observed"))$Observed
    bac.PA@sam_data$even <- bac.PA@sam_data$shannon/log(bac.PA@sam_data$richness)
    PAsample.data.bac <- data.frame(bac.PA@sam_data)
    PAsample.data.bac$Time <- factor(PAsample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    PA_otu <- bac.PA@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    PAphylo.diversity <- pd(t(PA_otu), bac.PA@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    PAphylo.diversity$Samples <- rownames(PAphylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    PAphylo.diversity.join <- left_join(PAphylo.diversity, as.data.frame(bac.PA@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    PAphylo.diversity.join$Time.Point <- factor(PAphylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## Arkansas Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.AR@sam_data$shannon <- estimate_richness(bac.AR, measures=c("Shannon"))$Shannon
    bac.AR@sam_data$invsimpson <- estimate_richness(bac.AR, measures=c("InvSimpson"))$InvSimpson
    bac.AR@sam_data$simpson <- estimate_richness(bac.AR, measures=c("Simpson"))$Simpson
    bac.AR@sam_data$richness <- estimate_richness(bac.AR, measures=c("Observed"))$Observed
    bac.AR@sam_data$even <- bac.AR@sam_data$shannon/log(bac.AR@sam_data$richness)
    ARsample.data.bac <- data.frame(bac.AR@sam_data)
    ARsample.data.bac$Time <- factor(ARsample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    AR_otu <- bac.AR@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    ARphylo.diversity <- pd(t(AR_otu), bac.AR@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    ARphylo.diversity$Samples <- rownames(ARphylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    ARphylo.diversity.join <- left_join(ARphylo.diversity, as.data.frame(bac.AR@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    ARphylo.diversity.join$Time.Point <- factor(ARphylo.diversity.join$Time, levels = c("0", "8", "12","16"))

## TVREC Soil

    ### Calculating Shannon, Simpson, Inverse Simpson and Pielou's Eveness ###

    bac.TVR@sam_data$shannon <- estimate_richness(bac.TVR, measures=c("Shannon"))$Shannon
    bac.TVR@sam_data$invsimpson <- estimate_richness(bac.TVR, measures=c("InvSimpson"))$InvSimpson
    bac.TVR@sam_data$simpson <- estimate_richness(bac.TVR, measures=c("Simpson"))$Simpson
    bac.TVR@sam_data$richness <- estimate_richness(bac.TVR, measures=c("Observed"))$Observed
    bac.TVR@sam_data$even <- bac.TVR@sam_data$shannon/log(bac.TVR@sam_data$richness)
    TVRsample.data.bac <- data.frame(bac.TVR@sam_data)
    TVRsample.data.bac$Time <- factor(TVRsample.data.bac$Time, levels = c("0", "8", "12","16"))
    ### Calculating Faiths Phylogenetic Diversity ###
    TVR_otu <- bac.TVR@otu_table %>%
      as.data.frame()

    # Calculate Faith's phylogenetic diversity for each sample, excluding the root
    TVRphylo.diversity <- pd(t(TVR_otu), bac.TVR@phy_tree, include.root = F)

    # Add sample names as a column in the phylogenetic diversity data frame
    TVRphylo.diversity$Samples <- rownames(TVRphylo.diversity)

    # Join the phylogenetic diversity data frame with the sample data from 'bac.Nonsteamed' based on sample names
    TVRphylo.diversity.join <- left_join(TVRphylo.diversity, as.data.frame(bac.TVR@sam_data), by = c("Samples" = "Sample_ID"))

    # Convert the 'Time' column in the joined data frame to a factor with specified levels
    TVRphylo.diversity.join$Time.Point <- factor(TVRphylo.diversity.join$Time, levels = c("0", "8", "12","16"))

\###Merging Steamed and Nonsteamed Alpha Diversity Data Frames###

    # Creating new columns To distinguish the dataframes
    autophylo.diversity.join$Location2 <- "Autoclaved_Control"
    MIphylo.diversity.join$Location2 <- "Michigan"
    NDphylo.diversity.join$Location2 <- "North_Dakota"
    Wirephylo.diversity.join$Location2 <- "Wiregrass"
    PAphylo.diversity.join$Location2 <- "Pennsylvania"
    ARphylo.diversity.join$Location2 <- "Arkansas"
    EVphylo.diversity.join$Location2 <- "EV_Smith"
    TVRphylo.diversity.join$Location2 <- "TVREC"
    # Combine the 'Sbphylo.diversity.join' and 'nbphylo.diversity.join' data frames into one data frame 'BAlphaCombined'
    BAlphaCombined <- rbind.data.frame(autophylo.diversity.join, MIphylo.diversity.join,NDphylo.diversity.join, Wirephylo.diversity.join, PAphylo.diversity.join, ARphylo.diversity.join, EVphylo.diversity.join, TVRphylo.diversity.join)

\###Plotting Richness for Prokaryotes as a Box-Plot###

    bac.richness <- BAlphaCombined %>%
      ggplot(aes(x = Sample_Type, y = richness, color = Time)) +
      facet_wrap(~factor(Location2)) +
      geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
      geom_point(position=position_jitterdodge(0.05)) +
      geom_pwc(aes(group = Location2), method = "t_test", label = "p.adj.format") +
      #stat_compare_means(method = "t.test") +
      scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)) +
      scale_fill_manual(values = cbbPalette ) +
      color_palette(cbbPalette) +
      #stat_summary(fun = mean,geom="line") +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
      ylab("Richness") +
      xlab("Type of Environment")+
      #geom_jitter(width = 0.1, alpha = 0.8)+
      #scale_color_manual(values=cbbPalette, name = "", labels = c("Planting", "17hrs", "V2")) +  
      #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
      theme_classic()+
      theme(legend.text = element_text(face = "italic", size = 12), legend.title = element_blank())
    bac.richness

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Plotting%20Richness-1.png)

    #The code begins by loading the necessary libraries, dplyr, ggplot2, and ggpubr, which are essential for data manipulation and visualization. The bac.richness object is created using the BAlphaCombined data frame. The ggplot function is initialized with aesthetic mappings, where Sample_Type is mapped to the x-axis, richness to the y-axis, and Time to the color of the points. The facet_wrap(~factor(Location2)) function is used to create separate plots for each level of Location2, allowing for a clear comparison across different locations.

    #Boxplots are added to the plot using geom_boxplot(position = position_dodge2(0.85, preserve = "single")), which ensures that the boxplots are dodged to avoid overlap. To show individual data points, geom_point(position = position_jitterdodge(0.05)) is used, adding jittered points to the plot. Pairwise comparisons are added using geom_pwc(aes(group = Location2), method = "t_test", label = "p.adj.format"), which performs t-tests and adjusts the p-values for multiple comparisons.

    #The y-axis limits and breaks are set using scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)), ensuring the y-axis ranges from 0 to 5000 with breaks at every 500 units. The fill colors for the plot are manually set using scale_fill_manual(values = cbbPalette), and the color palette for points and lines is defined with color_palette(cbbPalette). The y-axis and x-axis are labeled "Richness" and "Type of Environment" respectively, using ylab("Richness") and xlab("Type of Environment").

    #A classic theme is applied to the plot with theme_classic(), giving it a clean and professional look. Finally, the legend text is customized to be italic and sized at 12, and the legend title is removed using theme(legend.text = element_text(face = "italic", size = 12), legend.title = element_blank()). The plot is then printed by calling bac.richness, displaying the final visualization.

    bac.richness <- BAlphaCombined %>%
      ggplot(aes(x = Sample_Type, y = richness, color = Location2)) +
      facet_wrap(~factor(Time, levels = c("0", "8", "12","16"))) +
      geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
      geom_point(position=position_jitterdodge(0.05)) +
      geom_pwc(aes(group = Location2), method = "t_test", label = "p.adj.format") +
      #stat_compare_means(method = "t.test") +
      scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)) +
      scale_fill_manual(values = cbbPalette ) +
      color_palette(cbbPalette) +
      #stat_summary(fun = mean,geom="line") +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
      ylab("Richness") +
      xlab("Type of Environment")+
      #geom_jitter(width = 0.1, alpha = 0.8)+
      #scale_color_manual(values=cbbPalette, name = "", labels = c("Planting", "17hrs", "V2")) +  
      #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
      theme_classic()+
      theme(legend.text = element_text(face = "italic", size = 12), legend.title = element_blank())
    bac.richness

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Plotting%20Richness-2.png)

    #The code begins by loading the necessary libraries, dplyr, ggplot2, and ggpubr, which are essential for data manipulation and visualization. The bac.richness object is created using the BAlphaCombined data frame. The ggplot function is initialized with aesthetic mappings, where Sample_Type is mapped to the x-axis, richness to the y-axis, and Location2 to the color of the points. The facet_wrap(~factor(Time, levels = c("0", "8", "12","16"))) function is used to create separate plots for each level of Time, allowing for a clear comparison across different time points.

    #Boxplots are added to the plot using geom_boxplot(position = position_dodge2(0.85, preserve = "single")), which ensures that the boxplots are dodged to avoid overlap. To show individual data points, geom_point(position = position_jitterdodge(0.05)) is used, adding jittered points to the plot. Pairwise comparisons are added using geom_pwc(aes(group = Location2), method = "t_test", label = "p.adj.format"), which performs t-tests and adjusts the p-values for multiple comparisons.

    #The y-axis limits and breaks are set using scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)), ensuring the y-axis ranges from 0 to 5000 with breaks at every 500 units. The fill colors for the plot are manually set using scale_fill_manual(values = cbbPalette), and the color palette for points and lines is defined with color_palette(cbbPalette). The y-axis and x-axis are labeled "Richness" and "Type of Environment" respectively, using ylab("Richness") and xlab("Type of Environment").

    #A classic theme is applied to the plot with theme_classic(), giving it a clean and professional look. Finally, the legend text is customized to be italic and sized at 12, and the legend title is removed using theme(legend.text = element_text(face = "italic", size = 12), legend.title = element_blank()). The plot is then printed by calling bac.richness, displaying the final visualization.

\###Plotting Faith’s Phylogenetic Diversity for Prokaryotes as a
Box-Plot###

    bac.phylo.div <- BAlphaCombined %>%
      ggplot(aes(x = Sample_Type, y = PD, color = Location2)) +
      facet_wrap(~factor(Time, levels = c("0", "8", "12","16")), scale = "free") +
      #stat_compare_means(method = "t.test") +
      geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
        #scale_y_continuous(limits = c(30, 210), breaks = seq(30, 210, 70)) +
      geom_point(position=position_jitterdodge(0.05)) +
      scale_fill_manual(values = cbbPalette ) +
      color_palette(cbbPalette) +
      #stat_summary(fun = mean,geom="line") +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
      ylab("Faith's PD") +
      xlab("Type of Environment")+
      #geom_jitter(width = 0.1, alpha = 0.8)+
      #scale_color_manual(values=cbbPalette, name = "", labels = c("Planting", "17hrs", "V2")) +  
      #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
      theme_classic()+
      theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank())
    bac.phylo.div

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Plotting%20Faiths%20Phylogenetic%20Diversity-1.png)

    #The code begins by loading the necessary libraries, dplyr, ggplot2, and ggpubr, which are essential for data manipulation and visualization. The bac.phylo.div object is created using the BAlphaCombined data frame. The ggplot function is initialized with aesthetic mappings, where Sample_Type is mapped to the x-axis, PD to the y-axis, and Location2 to the color of the points. The facet_wrap(~factor(Time, levels = c("0", "8", "12","16")), scale = "free") function is used to create separate plots for each level of Time, allowing for a clear comparison across different time points, with each facet having its own scale.

    #Boxplots are added to the plot using geom_boxplot(position = position_dodge2(0.85, preserve = "single")), which ensures that the boxplots are dodged to avoid overlap. To show individual data points, geom_point(position = position_jitterdodge(0.05)) is used, adding jittered points to the plot. The y-axis limits and breaks are commented out, but if needed, they can be set using scale_y_continuous(limits = c(30, 210), breaks = seq(30, 210, 70)).

    #The fill colors for the plot are manually set using scale_fill_manual(values = cbbPalette), and the color palette for points and lines is defined with color_palette(cbbPalette). The y-axis and x-axis are labeled "Faith's PD" and "Type of Environment" respectively, using ylab("Faith's PD") and xlab("Type of Environment").

    #A classic theme is applied to the plot with theme_classic(), giving it a clean and professional look. Finally, the legend text is customized to be italic and sized at 8, and the legend title is removed using theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank()). The plot is then printed by calling bac.phylo.div, displaying the final visualization.

    bac.phylo.div <- BAlphaCombined %>%
      ggplot(aes(x = Sample_Type, y = PD, color = Time)) +
      facet_wrap(~factor(Location2), scale = "free") +
      #stat_compare_means(method = "t.test") +
      geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
        #scale_y_continuous(limits = c(30, 210), breaks = seq(30, 210, 70)) +
      geom_point(position=position_jitterdodge(0.05)) +
      scale_fill_manual(values = cbbPalette ) +
      color_palette(cbbPalette) +
      #stat_summary(fun = mean,geom="line") +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
      ylab("Faith's PD") +
      xlab("Type of Environment")+
      #geom_jitter(width = 0.1, alpha = 0.8)+
      #scale_color_manual(values=cbbPalette, name = "", labels = c("Planting", "17hrs", "V2")) +  
      #stat_compare_means(method = "kruskal", hide.ns = TRUE) +
      theme_classic()+
      theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank())
    bac.phylo.div

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Plotting%20Faiths%20Phylogenetic%20Diversity-2.png)

    #The code begins by loading the necessary libraries, dplyr, ggplot2, and ggpubr, which are essential for data manipulation and visualization. The bac.phylo.div object is created using the BAlphaCombined data frame. The ggplot function is initialized with aesthetic mappings, where Sample_Type is mapped to the x-axis, PD to the y-axis, and Time to the color of the points. The facet_wrap(~factor(Location2), scale = "free") function is used to create separate plots for each level of Location2, allowing for a clear comparison across different locations, with each facet having its own scale.

    #Boxplots are added to the plot using geom_boxplot(position = position_dodge2(0.85, preserve = "single")), which ensures that the boxplots are dodged to avoid overlap. To show individual data points, geom_point(position = position_jitterdodge(0.05)) is used, adding jittered points to the plot. The y-axis limits and breaks are commented out, but if needed, they can be set using scale_y_continuous(limits = c(30, 210), breaks = seq(30, 210, 70)).

    #The fill colors for the plot are manually set using scale_fill_manual(values = cbbPalette), and the color palette for points and lines is defined with color_palette(cbbPalette). The y-axis and x-axis are labeled "Faith's PD" and "Type of Environment" respectively, using ylab("Faith's PD") and xlab("Type of Environment").

    #A classic theme is applied to the plot with theme_classic(), giving it a clean and professional look. Finally, the legend text is customized to be italic and sized at 8, and the legend title is removed using theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank()). The plot is then printed by calling bac.phylo.div, displaying the final visualization.

\###Plotting Pielou’s Eveness for Prokaryotes as a Box-Plot###

    bac.even <- BAlphaCombined %>%
      ggplot(aes(x = Sample_Type, y = even, color = Location2)) +
      facet_wrap(~factor(Time, levels = c("0", "8", "12","16")), scale = "free") +
      #scale_y_continuous(limits = c(0.45, 0.85), breaks = seq(0.45, 0.85, 0.05)) +
      #stat_compare_means(method = "t.test") +
      geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
      geom_point(position=position_jitterdodge(0.05)) +
      scale_fill_manual(values = cbbPalette ) +
      color_palette(cbbPalette) +
      #stat_summary(fun.y=mean,geom="bar", position = "dodge", width = 0.5) +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
      ylab("Pielou's evenness") +
      xlab("Type of Environment")+
      #scale_color_manual(values=cbbPalette, name = "", labels = c("Planting", "17hrs", "V2")) +  
      theme_classic() +
      theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank())
    bac.even

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Plotting%20Pielous%20Eveness-1.png)

    #The code begins by loading the necessary libraries, dplyr, ggplot2, and ggpubr, which are essential for data manipulation and visualization. The bac.even object is created using the BAlphaCombined data frame. The ggplot function is initialized with aesthetic mappings, where Sample_Type is mapped to the x-axis, even to the y-axis, and Location2 to the color of the points. The facet_wrap(~factor(Time, levels = c("0", "8", "12","16")), scale = "free") function is used to create separate plots for each level of Time, allowing for a clear comparison across different time points, with each facet having its own scale.

    #Boxplots are added to the plot using geom_boxplot(position = position_dodge2(0.85, preserve = "single")), which ensures that the boxplots are dodged to avoid overlap. To show individual data points, geom_point(position = position_jitterdodge(0.05)) is used, adding jittered points to the plot. The y-axis limits and breaks are commented out, but if needed, they can be set using scale_y_continuous(limits = c(0.45, 0.85), breaks = seq(0.45, 0.85, 0.05)).

    #The fill colors for the plot are manually set using scale_fill_manual(values = cbbPalette), and the color palette for points and lines is defined with color_palette(cbbPalette). The y-axis and x-axis are labeled "Pielou's evenness" and "Type of Environment" respectively, using ylab("Pielou's evenness") and xlab("Type of Environment").

    #A classic theme is applied to the plot with theme_classic(), giving it a clean and professional look. Finally, the legend text is customized to be italic and sized at 8, and the legend title is removed using theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank()). The plot is then printed by calling bac.even, displaying the final visualization.

    bac.even <- BAlphaCombined %>%
      ggplot(aes(x = Sample_Type, y = even, color = Time)) +
      facet_wrap(~factor(Location2)) +
      #scale_y_continuous(limits = c(0.45, 0.85), breaks = seq(0.45, 0.85, 0.05)) +
      #stat_compare_means(method = "t.test") +
      geom_boxplot(position = position_dodge2(0.85, preserve = "single")) + 
      geom_point(position=position_jitterdodge(0.05)) +
      scale_fill_manual(values = cbbPalette ) +
      color_palette(cbbPalette) +
      #stat_summary(fun.y=mean,geom="bar", position = "dodge", width = 0.5) +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
      ylab("Pielou's evenness") +
      xlab("Type of Environment")+
      #scale_color_manual(values=cbbPalette, name = "", labels = c("Planting", "17hrs", "V2")) +  
      theme_classic() +
      theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank())
    bac.even

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Plotting%20Pielous%20Eveness-2.png)

    #The code begins by loading the necessary libraries, dplyr, ggplot2, and ggpubr, which are essential for data manipulation and visualization. The bac.even object is created using the BAlphaCombined data frame. The ggplot function is initialized with aesthetic mappings, where Sample_Type is mapped to the x-axis, even to the y-axis, and Time to the color of the points. The facet_wrap(~factor(Location2)) function is used to create separate plots for each level of Location2, allowing for a clear comparison across different locations.

    #Boxplots are added to the plot using geom_boxplot(position = position_dodge2(0.85, preserve = "single")), which ensures that the boxplots are dodged to avoid overlap. To show individual data points, geom_point(position = position_jitterdodge(0.05)) is used, adding jittered points to the plot. The y-axis limits and breaks are commented out, but if needed, they can be set using scale_y_continuous(limits = c(0.45, 0.85), breaks = seq(0.45, 0.85, 0.05)).

    #The fill colors for the plot are manually set using scale_fill_manual(values = cbbPalette), and the color palette for points and lines is defined with color_palette(cbbPalette). The y-axis and x-axis are labeled "Pielou's evenness" and "Type of Environment" respectively, using ylab("Pielou's evenness") and xlab("Type of Environment").

    #A classic theme is applied to the plot with theme_classic(), giving it a clean and professional look. Finally, the legend text is customized to be italic and sized at 8, and the legend title is removed using theme(legend.text = element_text(face = "italic", size = 8), legend.title = element_blank()). The plot is then printed by calling bac.even, displaying the final visualization.

#### Richness statistics

    # for every X - unit change in X we observed a exp^b times as many species.  
    # 1 exp^b
    # 2. Times as many for counts 
    # 3. report confidence interval - exp^lower to exp^upper

    ### calculating observed means ###

    ObservedRichness <- BAlphaCombined %>%
      group_by(Location_EC, Time, Sample_Type) %>%  
      summarise(MeanRichness = mean(richness), STDRichness = sd(richness)) %>% 
      arrange(-MeanRichness)

    ## `summarise()` has grouped output by 'Location_EC', 'Time'. You can override
    ## using the `.groups` argument.

    ### Richness Statistics ###

    # Since the data is may not normally distributed because and it is count data. Fit a Poisson regression model to predict 'richness' based on 'Type.revised', 'Time', 'Trt', and their interactions
    RichPoisson <- glm(richness ~ Sample_Type + Time + Location2 + Location2 * Time * Sample_Type, data = BAlphaCombined, family = "poisson")

    # Perform an analysis of variance (ANOVA) on the fitted model
    anova(RichPoisson)

    ## Analysis of Deviance Table
    ## 
    ## Model: poisson, link: log
    ## 
    ## Response: richness
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##                            Df Deviance Resid. Df Resid. Dev
    ## NULL                                         383     161275
    ## Sample_Type                 1        5       382     161270
    ## Time                        3      160       379     161110
    ## Location2                   7   134126       372      26985
    ## Time:Location2             21     2449       351      24536
    ## Sample_Type:Location2       7     9212       344      15325
    ## Sample_Type:Time            3       51       341      15273
    ## Sample_Type:Time:Location2 21     3050       320      12223
    ##                                         Pr(>Chi)    
    ## NULL                                                
    ## Sample_Type                              0.02261 *  
    ## Time                       < 0.00000000000000022 ***
    ## Location2                  < 0.00000000000000022 ***
    ## Time:Location2             < 0.00000000000000022 ***
    ## Sample_Type:Location2      < 0.00000000000000022 ***
    ## Sample_Type:Time                 0.0000000000416 ***
    ## Sample_Type:Time:Location2 < 0.00000000000000022 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # Display a summary of the fitted model, including coefficients, standard errors, z-values, and p-values
    summary(RichPoisson)

    ## 
    ## Call:
    ## glm(formula = richness ~ Sample_Type + Time + Location2 + Location2 * 
    ##     Time * Sample_Type, family = "poisson", data = BAlphaCombined)
    ## 
    ## Coefficients:
    ##                                                              Estimate
    ## (Intercept)                                                 8.1364207
    ## Sample_TypeSpermosphere                                    -0.0504208
    ## Time12                                                      0.0377656
    ## Time16                                                      0.0735224
    ## Time8                                                      -0.0856110
    ## Location2Autoclaved_Control                                -1.3777456
    ## Location2EV_Smith                                           0.2018053
    ## Location2Michigan                                           0.1312422
    ## Location2North_Dakota                                      -0.3948871
    ## Location2Pennsylvania                                       0.0769159
    ## Location2TVREC                                              0.0506004
    ## Location2Wiregrass                                          0.1748957
    ## Time12:Location2Autoclaved_Control                         -0.1386436
    ## Time16:Location2Autoclaved_Control                         -0.4657270
    ## Time8:Location2Autoclaved_Control                          -0.2438823
    ## Time12:Location2EV_Smith                                    0.0004825
    ## Time16:Location2EV_Smith                                   -0.0613578
    ## Time8:Location2EV_Smith                                     0.1132125
    ## Time12:Location2Michigan                                   -0.0512031
    ## Time16:Location2Michigan                                   -0.0991310
    ## Time8:Location2Michigan                                     0.0814527
    ## Time12:Location2North_Dakota                               -0.0667063
    ## Time16:Location2North_Dakota                               -0.0997838
    ## Time8:Location2North_Dakota                                -0.0125195
    ## Time12:Location2Pennsylvania                               -0.0697548
    ## Time16:Location2Pennsylvania                               -0.0820505
    ## Time8:Location2Pennsylvania                                 0.0539949
    ## Time12:Location2TVREC                                       0.0324261
    ## Time16:Location2TVREC                                      -0.1004063
    ## Time8:Location2TVREC                                        0.1745181
    ## Time12:Location2Wiregrass                                  -0.0160921
    ## Time16:Location2Wiregrass                                  -0.0430710
    ## Time8:Location2Wiregrass                                    0.0773040
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control         0.0051068
    ## Sample_TypeSpermosphere:Location2EV_Smith                   0.0649849
    ## Sample_TypeSpermosphere:Location2Michigan                   0.0511905
    ## Sample_TypeSpermosphere:Location2North_Dakota               0.0339228
    ## Sample_TypeSpermosphere:Location2Pennsylvania               0.0608902
    ## Sample_TypeSpermosphere:Location2TVREC                      0.0605713
    ## Sample_TypeSpermosphere:Location2Wiregrass                  0.0806735
    ## Sample_TypeSpermosphere:Time12                             -0.1069527
    ## Sample_TypeSpermosphere:Time16                             -0.1581115
    ## Sample_TypeSpermosphere:Time8                              -0.0290730
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control  1.1396281
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control  1.2314481
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   1.1621824
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith            0.0746979
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith            0.1320603
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith            -0.0526547
    ## Sample_TypeSpermosphere:Time12:Location2Michigan            0.0037714
    ## Sample_TypeSpermosphere:Time16:Location2Michigan            0.1163533
    ## Sample_TypeSpermosphere:Time8:Location2Michigan            -0.0009459
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota        0.1958164
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota        0.1309311
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota         0.1504122
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania        0.1141878
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania        0.0244117
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania         0.0669274
    ## Sample_TypeSpermosphere:Time12:Location2TVREC               0.0215430
    ## Sample_TypeSpermosphere:Time16:Location2TVREC               0.1676777
    ## Sample_TypeSpermosphere:Time8:Location2TVREC               -0.1470792
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass           0.0650506
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass           0.0955217
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass            0.0450980
    ##                                                            Std. Error  z value
    ## (Intercept)                                                 0.0069843 1164.958
    ## Sample_TypeSpermosphere                                     0.0100042   -5.040
    ## Time12                                                      0.0097854    3.859
    ## Time16                                                      0.0097007    7.579
    ## Time8                                                       0.0100956   -8.480
    ## Location2Autoclaved_Control                                 0.0155641  -88.521
    ## Location2EV_Smith                                           0.0094152   21.434
    ## Location2Michigan                                           0.0095688   13.716
    ## Location2North_Dakota                                       0.0110082  -35.872
    ## Location2Pennsylvania                                       0.0096928    7.935
    ## Location2TVREC                                              0.0097547    5.187
    ## Location2Wiregrass                                          0.0094728   18.463
    ## Time12:Location2Autoclaved_Control                          0.0224323   -6.181
    ## Time16:Location2Autoclaved_Control                          0.0239569  -19.440
    ## Time8:Location2Autoclaved_Control                           0.0237559  -10.266
    ## Time12:Location2EV_Smith                                    0.0131905    0.037
    ## Time16:Location2EV_Smith                                    0.0131664   -4.660
    ## Time8:Location2EV_Smith                                     0.0134376    8.425
    ## Time12:Location2Michigan                                    0.0134868   -3.797
    ## Time16:Location2Michigan                                    0.0134452   -7.373
    ## Time8:Location2Michigan                                     0.0136989    5.946
    ## Time12:Location2North_Dakota                                0.0155782   -4.282
    ## Time16:Location2North_Dakota                                0.0155187   -6.430
    ## Time8:Location2North_Dakota                                 0.0159433   -0.785
    ## Time12:Location2Pennsylvania                                0.0136952   -5.093
    ## Time16:Location2Pennsylvania                                0.0135952   -6.035
    ## Time8:Location2Pennsylvania                                 0.0139180    3.880
    ## Time12:Location2TVREC                                       0.0136146    2.382
    ## Time16:Location2TVREC                                       0.0137155   -7.321
    ## Time8:Location2TVREC                                        0.0138103   12.637
    ## Time12:Location2Wiregrass                                   0.0132959   -1.210
    ## Time16:Location2Wiregrass                                   0.0132205   -3.258
    ## Time8:Location2Wiregrass                                    0.0135709    5.696
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control         0.0222704    0.229
    ## Sample_TypeSpermosphere:Location2EV_Smith                   0.0133880    4.854
    ## Sample_TypeSpermosphere:Location2Michigan                   0.0136240    3.757
    ## Sample_TypeSpermosphere:Location2North_Dakota               0.0156872    2.162
    ## Sample_TypeSpermosphere:Location2Pennsylvania               0.0137823    4.418
    ## Sample_TypeSpermosphere:Location2TVREC                      0.0138695    4.367
    ## Sample_TypeSpermosphere:Location2Wiregrass                  0.0134451    6.000
    ## Sample_TypeSpermosphere:Time12                              0.0142139   -7.525
    ## Sample_TypeSpermosphere:Time16                              0.0141858  -11.146
    ## Sample_TypeSpermosphere:Time8                               0.0145176   -2.003
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control  0.0298603   38.165
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control  0.0313995   39.219
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   0.0310844   37.388
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith            0.0189371    3.945
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith            0.0189636    6.964
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith             0.0192399   -2.737
    ## Sample_TypeSpermosphere:Time12:Location2Michigan            0.0194668    0.194
    ## Sample_TypeSpermosphere:Time16:Location2Michigan            0.0194006    5.997
    ## Sample_TypeSpermosphere:Time8:Location2Michigan             0.0195827   -0.048
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota        0.0221786    8.829
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota        0.0223432    5.860
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota         0.0225511    6.670
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania        0.0196051    5.824
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania        0.0197026    1.239
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania         0.0197907    3.382
    ## Sample_TypeSpermosphere:Time12:Location2TVREC               0.0195998    1.099
    ## Sample_TypeSpermosphere:Time16:Location2TVREC               0.0196942    8.514
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                0.0198887   -7.395
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass           0.0190615    3.413
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass           0.0190442    5.016
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass            0.0192917    2.338
    ##                                                                        Pr(>|z|)
    ## (Intercept)                                                < 0.0000000000000002
    ## Sample_TypeSpermosphere                                      0.0000004656043676
    ## Time12                                                                 0.000114
    ## Time16                                                       0.0000000000000348
    ## Time8                                                      < 0.0000000000000002
    ## Location2Autoclaved_Control                                < 0.0000000000000002
    ## Location2EV_Smith                                          < 0.0000000000000002
    ## Location2Michigan                                          < 0.0000000000000002
    ## Location2North_Dakota                                      < 0.0000000000000002
    ## Location2Pennsylvania                                        0.0000000000000021
    ## Location2TVREC                                               0.0000002133833774
    ## Location2Wiregrass                                         < 0.0000000000000002
    ## Time12:Location2Autoclaved_Control                           0.0000000006388095
    ## Time16:Location2Autoclaved_Control                         < 0.0000000000000002
    ## Time8:Location2Autoclaved_Control                          < 0.0000000000000002
    ## Time12:Location2EV_Smith                                               0.970820
    ## Time16:Location2EV_Smith                                     0.0000031590987378
    ## Time8:Location2EV_Smith                                    < 0.0000000000000002
    ## Time12:Location2Michigan                                               0.000147
    ## Time16:Location2Michigan                                     0.0000000000001669
    ## Time8:Location2Michigan                                      0.0000000027492311
    ## Time12:Location2North_Dakota                                 0.0000185199035371
    ## Time16:Location2North_Dakota                                 0.0000000001276988
    ## Time8:Location2North_Dakota                                            0.432305
    ## Time12:Location2Pennsylvania                                 0.0000003517785034
    ## Time16:Location2Pennsylvania                                 0.0000000015869417
    ## Time8:Location2Pennsylvania                                            0.000105
    ## Time12:Location2TVREC                                                  0.017232
    ## Time16:Location2TVREC                                        0.0000000000002468
    ## Time8:Location2TVREC                                       < 0.0000000000000002
    ## Time12:Location2Wiregrass                                              0.226162
    ## Time16:Location2Wiregrass                                              0.001122
    ## Time8:Location2Wiregrass                                     0.0000000122431180
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control                    0.818627
    ## Sample_TypeSpermosphere:Location2EV_Smith                    0.0000012102271167
    ## Sample_TypeSpermosphere:Location2Michigan                              0.000172
    ## Sample_TypeSpermosphere:Location2North_Dakota                          0.030584
    ## Sample_TypeSpermosphere:Location2Pennsylvania                0.0000099617467965
    ## Sample_TypeSpermosphere:Location2TVREC                       0.0000125829794466
    ## Sample_TypeSpermosphere:Location2Wiregrass                   0.0000000019703692
    ## Sample_TypeSpermosphere:Time12                               0.0000000000000529
    ## Sample_TypeSpermosphere:Time16                             < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time8                                          0.045220
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control  < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith             0.0000799596407376
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith             0.0000000000033099
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith                        0.006205
    ## Sample_TypeSpermosphere:Time12:Location2Michigan                       0.846385
    ## Sample_TypeSpermosphere:Time16:Location2Michigan             0.0000000020049566
    ## Sample_TypeSpermosphere:Time8:Location2Michigan                        0.961475
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota       < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota         0.0000000046285308
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota          0.0000000000256064
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania         0.0000000057317789
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania                   0.215342
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania                    0.000720
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                          0.271705
    ## Sample_TypeSpermosphere:Time16:Location2TVREC              < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                 0.0000000000001413
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass                      0.000643
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass            0.0000005281373170
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass                       0.019403
    ##                                                               
    ## (Intercept)                                                ***
    ## Sample_TypeSpermosphere                                    ***
    ## Time12                                                     ***
    ## Time16                                                     ***
    ## Time8                                                      ***
    ## Location2Autoclaved_Control                                ***
    ## Location2EV_Smith                                          ***
    ## Location2Michigan                                          ***
    ## Location2North_Dakota                                      ***
    ## Location2Pennsylvania                                      ***
    ## Location2TVREC                                             ***
    ## Location2Wiregrass                                         ***
    ## Time12:Location2Autoclaved_Control                         ***
    ## Time16:Location2Autoclaved_Control                         ***
    ## Time8:Location2Autoclaved_Control                          ***
    ## Time12:Location2EV_Smith                                      
    ## Time16:Location2EV_Smith                                   ***
    ## Time8:Location2EV_Smith                                    ***
    ## Time12:Location2Michigan                                   ***
    ## Time16:Location2Michigan                                   ***
    ## Time8:Location2Michigan                                    ***
    ## Time12:Location2North_Dakota                               ***
    ## Time16:Location2North_Dakota                               ***
    ## Time8:Location2North_Dakota                                   
    ## Time12:Location2Pennsylvania                               ***
    ## Time16:Location2Pennsylvania                               ***
    ## Time8:Location2Pennsylvania                                ***
    ## Time12:Location2TVREC                                      *  
    ## Time16:Location2TVREC                                      ***
    ## Time8:Location2TVREC                                       ***
    ## Time12:Location2Wiregrass                                     
    ## Time16:Location2Wiregrass                                  ** 
    ## Time8:Location2Wiregrass                                   ***
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control           
    ## Sample_TypeSpermosphere:Location2EV_Smith                  ***
    ## Sample_TypeSpermosphere:Location2Michigan                  ***
    ## Sample_TypeSpermosphere:Location2North_Dakota              *  
    ## Sample_TypeSpermosphere:Location2Pennsylvania              ***
    ## Sample_TypeSpermosphere:Location2TVREC                     ***
    ## Sample_TypeSpermosphere:Location2Wiregrass                 ***
    ## Sample_TypeSpermosphere:Time12                             ***
    ## Sample_TypeSpermosphere:Time16                             ***
    ## Sample_TypeSpermosphere:Time8                              *  
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control  ***
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith           ***
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith           ***
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith            ** 
    ## Sample_TypeSpermosphere:Time12:Location2Michigan              
    ## Sample_TypeSpermosphere:Time16:Location2Michigan           ***
    ## Sample_TypeSpermosphere:Time8:Location2Michigan               
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota       ***
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota       ***
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota        ***
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania       ***
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania          
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania        ***
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time16:Location2TVREC              ***
    ## Sample_TypeSpermosphere:Time8:Location2TVREC               ***
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass          ***
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass          ***
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass           *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 161275  on 383  degrees of freedom
    ## Residual deviance:  12223  on 320  degrees of freedom
    ## AIC: 16125
    ## 
    ## Number of Fisher Scoring iterations: 4

    # Calculate the variance by dividing the deviance by the residual degrees of freedom
    RichPoisson$deviance / RichPoisson$df.residual

    ## [1] 38.19676

    #output is 38.20, way higher than what is recommended so negative binomial may be better, the target is to get close to 1.

    # To fix the distribution variance, Fit a negative binomial regression model to predict 'richness' based on 'Type.revised', 'Time', 'Trt', and their interactions
    RichStat <- glm.nb(richness ~ Sample_Type + Time + Location2 + Location2 * Time * Sample_Type, data = BAlphaCombined)
    # Perform an analysis of variance (ANOVA) on the fitted model
    anova(RichStat)

    ## Warning in anova.negbin(RichStat): tests made without re-estimating 'theta'

    ## Analysis of Deviance Table
    ## 
    ## Model: Negative Binomial(45.5464), link: log
    ## 
    ## Response: richness
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##                            Df Deviance Resid. Df Resid. Dev
    ## NULL                                         383     3240.0
    ## Sample_Type                 1     0.07       382     3239.9
    ## Time                        3     2.20       379     3237.7
    ## Location2                   7  2357.17       372      880.5
    ## Time:Location2             21    72.78       351      807.7
    ## Sample_Type:Location2       7   288.29       344      519.5
    ## Sample_Type:Time            3     7.70       341      511.8
    ## Sample_Type:Time:Location2 21   115.68       320      396.1
    ##                                         Pr(>Chi)    
    ## NULL                                                
    ## Sample_Type                              0.78868    
    ## Time                                     0.53158    
    ## Location2                  < 0.00000000000000022 ***
    ## Time:Location2              0.000000125195350095 ***
    ## Sample_Type:Location2      < 0.00000000000000022 ***
    ## Sample_Type:Time                         0.05255 .  
    ## Sample_Type:Time:Location2  0.000000000000004397 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # Display a summary of the fitted model, including coefficients, standard errors, z-values, and p-values
    summary(RichStat)

    ## 
    ## Call:
    ## glm.nb(formula = richness ~ Sample_Type + Time + Location2 + 
    ##     Location2 * Time * Sample_Type, data = BAlphaCombined, init.theta = 45.54644757, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                                              Estimate
    ## (Intercept)                                                 8.1364207
    ## Sample_TypeSpermosphere                                    -0.0504208
    ## Time12                                                      0.0377656
    ## Time16                                                      0.0735224
    ## Time8                                                      -0.0856110
    ## Location2Autoclaved_Control                                -1.3777456
    ## Location2EV_Smith                                           0.2018053
    ## Location2Michigan                                           0.1312422
    ## Location2North_Dakota                                      -0.3948871
    ## Location2Pennsylvania                                       0.0769159
    ## Location2TVREC                                              0.0506004
    ## Location2Wiregrass                                          0.1748957
    ## Time12:Location2Autoclaved_Control                         -0.1386436
    ## Time16:Location2Autoclaved_Control                         -0.4657270
    ## Time8:Location2Autoclaved_Control                          -0.2438823
    ## Time12:Location2EV_Smith                                    0.0004825
    ## Time16:Location2EV_Smith                                   -0.0613578
    ## Time8:Location2EV_Smith                                     0.1132125
    ## Time12:Location2Michigan                                   -0.0512031
    ## Time16:Location2Michigan                                   -0.0991310
    ## Time8:Location2Michigan                                     0.0814527
    ## Time12:Location2North_Dakota                               -0.0667063
    ## Time16:Location2North_Dakota                               -0.0997838
    ## Time8:Location2North_Dakota                                -0.0125195
    ## Time12:Location2Pennsylvania                               -0.0697548
    ## Time16:Location2Pennsylvania                               -0.0820505
    ## Time8:Location2Pennsylvania                                 0.0539949
    ## Time12:Location2TVREC                                       0.0324261
    ## Time16:Location2TVREC                                      -0.1004063
    ## Time8:Location2TVREC                                        0.1745181
    ## Time12:Location2Wiregrass                                  -0.0160921
    ## Time16:Location2Wiregrass                                  -0.0430710
    ## Time8:Location2Wiregrass                                    0.0773040
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control         0.0051068
    ## Sample_TypeSpermosphere:Location2EV_Smith                   0.0649849
    ## Sample_TypeSpermosphere:Location2Michigan                   0.0511905
    ## Sample_TypeSpermosphere:Location2North_Dakota               0.0339228
    ## Sample_TypeSpermosphere:Location2Pennsylvania               0.0608902
    ## Sample_TypeSpermosphere:Location2TVREC                      0.0605713
    ## Sample_TypeSpermosphere:Location2Wiregrass                  0.0806735
    ## Sample_TypeSpermosphere:Time12                             -0.1069527
    ## Sample_TypeSpermosphere:Time16                             -0.1581115
    ## Sample_TypeSpermosphere:Time8                              -0.0290730
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control  1.1396281
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control  1.2314481
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   1.1621824
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith            0.0746979
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith            0.1320603
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith            -0.0526547
    ## Sample_TypeSpermosphere:Time12:Location2Michigan            0.0037714
    ## Sample_TypeSpermosphere:Time16:Location2Michigan            0.1163533
    ## Sample_TypeSpermosphere:Time8:Location2Michigan            -0.0009459
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota        0.1958164
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota        0.1309311
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota         0.1504122
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania        0.1141878
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania        0.0244117
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania         0.0669274
    ## Sample_TypeSpermosphere:Time12:Location2TVREC               0.0215430
    ## Sample_TypeSpermosphere:Time16:Location2TVREC               0.1676777
    ## Sample_TypeSpermosphere:Time8:Location2TVREC               -0.1470792
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass           0.0650506
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass           0.0955217
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass            0.0450980
    ##                                                            Std. Error z value
    ## (Intercept)                                                 0.0608937 133.617
    ## Sample_TypeSpermosphere                                     0.0861314  -0.585
    ## Time12                                                      0.0861063   0.439
    ## Time16                                                      0.0860967   0.854
    ## Time8                                                       0.0861421  -0.994
    ## Location2Autoclaved_Control                                 0.0869527 -15.845
    ## Location2EV_Smith                                           0.0860650   2.345
    ## Location2Michigan                                           0.0860819   1.525
    ## Location2North_Dakota                                       0.0862538  -4.578
    ## Location2Pennsylvania                                       0.0860958   0.893
    ## Location2TVREC                                              0.0861028   0.588
    ## Location2Wiregrass                                          0.0860713   2.032
    ## Time12:Location2Autoclaved_Control                          0.1230459  -1.127
    ## Time16:Location2Autoclaved_Control                          0.1233329  -3.776
    ## Time8:Location2Autoclaved_Control                           0.1232940  -1.978
    ## Time12:Location2EV_Smith                                    0.1217007   0.004
    ## Time16:Location2EV_Smith                                    0.1216981  -0.504
    ## Time8:Location2EV_Smith                                     0.1217277   0.930
    ## Time12:Location2Michigan                                    0.1217332  -0.421
    ## Time16:Location2Michigan                                    0.1217286  -0.814
    ## Time8:Location2Michigan                                     0.1217569   0.669
    ## Time12:Location2North_Dakota                                0.1219826  -0.547
    ## Time16:Location2North_Dakota                                0.1219750  -0.818
    ## Time8:Location2North_Dakota                                 0.1220298  -0.103
    ## Time12:Location2Pennsylvania                                0.1217564  -0.573
    ## Time16:Location2Pennsylvania                                0.1217452  -0.674
    ## Time8:Location2Pennsylvania                                 0.1217817   0.443
    ## Time12:Location2TVREC                                       0.1217474   0.266
    ## Time16:Location2TVREC                                       0.1217587  -0.825
    ## Time8:Location2TVREC                                        0.1217694   1.433
    ## Time12:Location2Wiregrass                                   0.1217122  -0.132
    ## Time16:Location2Wiregrass                                   0.1217040  -0.354
    ## Time8:Location2Wiregrass                                    0.1217425   0.635
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control         0.1230164   0.042
    ## Sample_TypeSpermosphere:Location2EV_Smith                   0.1217223   0.534
    ## Sample_TypeSpermosphere:Location2Michigan                   0.1217484   0.420
    ## Sample_TypeSpermosphere:Location2North_Dakota               0.1219966   0.278
    ## Sample_TypeSpermosphere:Location2Pennsylvania               0.1217663   0.500
    ## Sample_TypeSpermosphere:Location2TVREC                      0.1217762   0.497
    ## Sample_TypeSpermosphere:Location2Wiregrass                  0.1217286   0.663
    ## Sample_TypeSpermosphere:Time12                              0.1218159  -0.878
    ## Sample_TypeSpermosphere:Time16                              0.1218126  -1.298
    ## Sample_TypeSpermosphere:Time8                               0.1218517  -0.239
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control  0.1736830   6.562
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control  0.1739542   7.079
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   0.1738976   6.683
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith            0.1721417   0.434
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith            0.1721446   0.767
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith             0.1721753  -0.306
    ## Sample_TypeSpermosphere:Time12:Location2Michigan            0.1722008   0.022
    ## Sample_TypeSpermosphere:Time16:Location2Michigan            0.1721933   0.676
    ## Sample_TypeSpermosphere:Time8:Location2Michigan             0.1722139  -0.005
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota        0.1725284   1.135
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota        0.1725496   0.759
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota         0.1725766   0.872
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania        0.1722164   0.663
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania        0.1722276   0.142
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania         0.1722377   0.389
    ## Sample_TypeSpermosphere:Time12:Location2TVREC               0.1722158   0.125
    ## Sample_TypeSpermosphere:Time16:Location2TVREC               0.1722266   0.974
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                0.1722490  -0.854
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass           0.1721554   0.378
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass           0.1721535   0.555
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass            0.1721811   0.262
    ##                                                                        Pr(>|z|)
    ## (Intercept)                                                < 0.0000000000000002
    ## Sample_TypeSpermosphere                                                0.558283
    ## Time12                                                                 0.660956
    ## Time16                                                                 0.393132
    ## Time8                                                                  0.320303
    ## Location2Autoclaved_Control                                < 0.0000000000000002
    ## Location2EV_Smith                                                      0.019037
    ## Location2Michigan                                                      0.127354
    ## Location2North_Dakota                                          0.00000468994993
    ## Location2Pennsylvania                                                  0.371656
    ## Location2TVREC                                                         0.556751
    ## Location2Wiregrass                                                     0.042155
    ## Time12:Location2Autoclaved_Control                                     0.259842
    ## Time16:Location2Autoclaved_Control                                     0.000159
    ## Time8:Location2Autoclaved_Control                                      0.047923
    ## Time12:Location2EV_Smith                                               0.996837
    ## Time16:Location2EV_Smith                                               0.614134
    ## Time8:Location2EV_Smith                                                0.352347
    ## Time12:Location2Michigan                                               0.674034
    ## Time16:Location2Michigan                                               0.415438
    ## Time8:Location2Michigan                                                0.503510
    ## Time12:Location2North_Dakota                                           0.584481
    ## Time16:Location2North_Dakota                                           0.413319
    ## Time8:Location2North_Dakota                                            0.918285
    ## Time12:Location2Pennsylvania                                           0.566710
    ## Time16:Location2Pennsylvania                                           0.500342
    ## Time8:Location2Pennsylvania                                            0.657495
    ## Time12:Location2TVREC                                                  0.789978
    ## Time16:Location2TVREC                                                  0.409580
    ## Time8:Location2TVREC                                                   0.151805
    ## Time12:Location2Wiregrass                                              0.894815
    ## Time16:Location2Wiregrass                                              0.723414
    ## Time8:Location2Wiregrass                                               0.525442
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control                    0.966887
    ## Sample_TypeSpermosphere:Location2EV_Smith                              0.593426
    ## Sample_TypeSpermosphere:Location2Michigan                              0.674148
    ## Sample_TypeSpermosphere:Location2North_Dakota                          0.780963
    ## Sample_TypeSpermosphere:Location2Pennsylvania                          0.617034
    ## Sample_TypeSpermosphere:Location2TVREC                                 0.618908
    ## Sample_TypeSpermosphere:Location2Wiregrass                             0.507502
    ## Sample_TypeSpermosphere:Time12                                         0.379951
    ## Sample_TypeSpermosphere:Time16                                         0.194291
    ## Sample_TypeSpermosphere:Time8                                          0.811421
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control     0.00000000005325
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control     0.00000000000145
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control      0.00000000002339
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith                       0.664337
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith                       0.442994
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith                        0.759741
    ## Sample_TypeSpermosphere:Time12:Location2Michigan                       0.982527
    ## Sample_TypeSpermosphere:Time16:Location2Michigan                       0.499223
    ## Sample_TypeSpermosphere:Time8:Location2Michigan                        0.995618
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota                   0.256383
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota                   0.447971
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota                    0.383444
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania                   0.507300
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania                   0.887285
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania                    0.697590
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                          0.900450
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                          0.330261
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                           0.393174
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass                      0.705535
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass                      0.578988
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass                       0.793381
    ##                                                               
    ## (Intercept)                                                ***
    ## Sample_TypeSpermosphere                                       
    ## Time12                                                        
    ## Time16                                                        
    ## Time8                                                         
    ## Location2Autoclaved_Control                                ***
    ## Location2EV_Smith                                          *  
    ## Location2Michigan                                             
    ## Location2North_Dakota                                      ***
    ## Location2Pennsylvania                                         
    ## Location2TVREC                                                
    ## Location2Wiregrass                                         *  
    ## Time12:Location2Autoclaved_Control                            
    ## Time16:Location2Autoclaved_Control                         ***
    ## Time8:Location2Autoclaved_Control                          *  
    ## Time12:Location2EV_Smith                                      
    ## Time16:Location2EV_Smith                                      
    ## Time8:Location2EV_Smith                                       
    ## Time12:Location2Michigan                                      
    ## Time16:Location2Michigan                                      
    ## Time8:Location2Michigan                                       
    ## Time12:Location2North_Dakota                                  
    ## Time16:Location2North_Dakota                                  
    ## Time8:Location2North_Dakota                                   
    ## Time12:Location2Pennsylvania                                  
    ## Time16:Location2Pennsylvania                                  
    ## Time8:Location2Pennsylvania                                   
    ## Time12:Location2TVREC                                         
    ## Time16:Location2TVREC                                         
    ## Time8:Location2TVREC                                          
    ## Time12:Location2Wiregrass                                     
    ## Time16:Location2Wiregrass                                     
    ## Time8:Location2Wiregrass                                      
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control           
    ## Sample_TypeSpermosphere:Location2EV_Smith                     
    ## Sample_TypeSpermosphere:Location2Michigan                     
    ## Sample_TypeSpermosphere:Location2North_Dakota                 
    ## Sample_TypeSpermosphere:Location2Pennsylvania                 
    ## Sample_TypeSpermosphere:Location2TVREC                        
    ## Sample_TypeSpermosphere:Location2Wiregrass                    
    ## Sample_TypeSpermosphere:Time12                                
    ## Sample_TypeSpermosphere:Time16                                
    ## Sample_TypeSpermosphere:Time8                                 
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control  ***
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith              
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith              
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith               
    ## Sample_TypeSpermosphere:Time12:Location2Michigan              
    ## Sample_TypeSpermosphere:Time16:Location2Michigan              
    ## Sample_TypeSpermosphere:Time8:Location2Michigan               
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota          
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota          
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota           
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania          
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania          
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania           
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                  
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass             
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass             
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass              
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(45.5464) family taken to be 1)
    ## 
    ##     Null deviance: 3239.97  on 383  degrees of freedom
    ## Residual deviance:  396.07  on 320  degrees of freedom
    ## AIC: 5911.4
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  45.55 
    ##           Std. Err.:  3.43 
    ## 
    ##  2 x log-likelihood:  -5781.401

    # Calculate the variance by dividing the deviance by the residual degrees of freedom
    RichStat$deviance / RichStat$df.residual

    ## [1] 1.237732

    # Mean separation
    # Calculate estimated marginal means (least-squares means) for the 'RichStat' model, examining the effect of treatment within combinations of type and time
    lsmeans <- emmeans(RichStat, ~Location2 | Sample_Type * Time)

    # Perform multiple comparisons with Tukey adjustment and generate compact letter display (CLD) for the estimated marginal means
    library(multcompView)

    ## Warning: package 'multcompView' was built under R version 4.4.3

    Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE, Letters = letters)

#### Evenness Statstics

    ##Evenness Statistics
    #install.packages("betareg")
    library(betareg)
    # Fit a beta regression model to eveness using specified predictors
    EvenStat.B <- betareg(even ~ Sample_Type + Time + Location2 + Location2 * Time * Sample_Type, data = BAlphaCombined)
    #works because the data ranges 0-1

    # Perform an ANOVA on the fitted model
    #anova(EvenStat.B)

    # Display a summary of the fitted model
    summary(EvenStat.B)

    ## 
    ## Call:
    ## betareg(formula = even ~ Sample_Type + Time + Location2 + Location2 * 
    ##     Time * Sample_Type, data = BAlphaCombined)
    ## 
    ## Quantile residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -5.8832 -0.1991  0.0545  0.3074  4.7879 
    ## 
    ## Coefficients (mean model with logit link):
    ##                                                            Estimate Std. Error
    ## (Intercept)                                                 1.47215    0.09667
    ## Sample_TypeSpermosphere                                     0.00092    0.13672
    ## Time12                                                     -0.09059    0.13489
    ## Time16                                                     -0.16939    0.13345
    ## Time8                                                       0.03045    0.13734
    ## Location2Autoclaved_Control                                -2.35371    0.12747
    ## Location2EV_Smith                                           0.31944    0.14438
    ## Location2Michigan                                           0.05502    0.13787
    ## Location2North_Dakota                                      -0.47696    0.12880
    ## Location2Pennsylvania                                       0.10717    0.13904
    ## Location2TVREC                                              0.01022    0.13691
    ## Location2Wiregrass                                          0.05608    0.13790
    ## Time12:Location2Autoclaved_Control                          0.53810    0.17637
    ## Time16:Location2Autoclaved_Control                          0.81263    0.17472
    ## Time8:Location2Autoclaved_Control                           0.15801    0.17946
    ## Time12:Location2EV_Smith                                    0.10722    0.20330
    ## Time16:Location2EV_Smith                                    0.10980    0.20088
    ## Time8:Location2EV_Smith                                    -0.02864    0.20465
    ## Time12:Location2Michigan                                    0.04401    0.19301
    ## Time16:Location2Michigan                                    0.23545    0.19378
    ## Time8:Location2Michigan                                    -0.03588    0.19535
    ## Time12:Location2North_Dakota                                0.07963    0.18070
    ## Time16:Location2North_Dakota                                0.12466    0.17932
    ## Time8:Location2North_Dakota                                -0.07455    0.18225
    ## Time12:Location2Pennsylvania                                0.12815    0.19601
    ## Time16:Location2Pennsylvania                                0.14105    0.19393
    ## Time8:Location2Pennsylvania                                -0.03990    0.19693
    ## Time12:Location2TVREC                                       0.08792    0.19232
    ## Time16:Location2TVREC                                       0.20357    0.19186
    ## Time8:Location2TVREC                                       -0.09985    0.19309
    ## Time12:Location2Wiregrass                                   0.13925    0.19453
    ## Time16:Location2Wiregrass                                   0.24212    0.19393
    ## Time8:Location2Wiregrass                                   -0.07452    0.19480
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control         0.48639    0.17763
    ## Sample_TypeSpermosphere:Location2EV_Smith                  -0.11556    0.20208
    ## Sample_TypeSpermosphere:Location2Michigan                   0.06247    0.19601
    ## Sample_TypeSpermosphere:Location2North_Dakota              -0.06619    0.18160
    ## Sample_TypeSpermosphere:Location2Pennsylvania              -0.16384    0.19421
    ## Sample_TypeSpermosphere:Location2TVREC                      0.05161    0.19444
    ## Sample_TypeSpermosphere:Location2Wiregrass                 -0.15689    0.19279
    ## Sample_TypeSpermosphere:Time12                             -0.68001    0.18441
    ## Sample_TypeSpermosphere:Time16                             -0.53116    0.18378
    ## Sample_TypeSpermosphere:Time8                              -0.57890    0.18772
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control  1.78093    0.24636
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control -1.48076    0.25246
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   1.61847    0.24761
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith            0.39873    0.27713
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith            0.33966    0.27577
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith             0.53298    0.28149
    ## Sample_TypeSpermosphere:Time12:Location2Michigan           -0.22109    0.26356
    ## Sample_TypeSpermosphere:Time16:Location2Michigan           -0.16296    0.26594
    ## Sample_TypeSpermosphere:Time8:Location2Michigan             0.24884    0.26989
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota        0.52917    0.24922
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota        0.46517    0.24880
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota         0.49344    0.25163
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania        0.48538    0.26745
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania        0.26388    0.26514
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania         0.52013    0.26996
    ## Sample_TypeSpermosphere:Time12:Location2TVREC               0.05816    0.26361
    ## Sample_TypeSpermosphere:Time16:Location2TVREC              -0.04623    0.26406
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                0.39419    0.26783
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass           0.44552    0.26534
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass           0.23530    0.26492
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass            0.33464    0.26586
    ##                                                            z value
    ## (Intercept)                                                 15.229
    ## Sample_TypeSpermosphere                                      0.007
    ## Time12                                                      -0.672
    ## Time16                                                      -1.269
    ## Time8                                                        0.222
    ## Location2Autoclaved_Control                                -18.465
    ## Location2EV_Smith                                            2.213
    ## Location2Michigan                                            0.399
    ## Location2North_Dakota                                       -3.703
    ## Location2Pennsylvania                                        0.771
    ## Location2TVREC                                               0.075
    ## Location2Wiregrass                                           0.407
    ## Time12:Location2Autoclaved_Control                           3.051
    ## Time16:Location2Autoclaved_Control                           4.651
    ## Time8:Location2Autoclaved_Control                            0.881
    ## Time12:Location2EV_Smith                                     0.527
    ## Time16:Location2EV_Smith                                     0.547
    ## Time8:Location2EV_Smith                                     -0.140
    ## Time12:Location2Michigan                                     0.228
    ## Time16:Location2Michigan                                     1.215
    ## Time8:Location2Michigan                                     -0.184
    ## Time12:Location2North_Dakota                                 0.441
    ## Time16:Location2North_Dakota                                 0.695
    ## Time8:Location2North_Dakota                                 -0.409
    ## Time12:Location2Pennsylvania                                 0.654
    ## Time16:Location2Pennsylvania                                 0.727
    ## Time8:Location2Pennsylvania                                 -0.203
    ## Time12:Location2TVREC                                        0.457
    ## Time16:Location2TVREC                                        1.061
    ## Time8:Location2TVREC                                        -0.517
    ## Time12:Location2Wiregrass                                    0.716
    ## Time16:Location2Wiregrass                                    1.249
    ## Time8:Location2Wiregrass                                    -0.383
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control          2.738
    ## Sample_TypeSpermosphere:Location2EV_Smith                   -0.572
    ## Sample_TypeSpermosphere:Location2Michigan                    0.319
    ## Sample_TypeSpermosphere:Location2North_Dakota               -0.364
    ## Sample_TypeSpermosphere:Location2Pennsylvania               -0.844
    ## Sample_TypeSpermosphere:Location2TVREC                       0.265
    ## Sample_TypeSpermosphere:Location2Wiregrass                  -0.814
    ## Sample_TypeSpermosphere:Time12                              -3.688
    ## Sample_TypeSpermosphere:Time16                              -2.890
    ## Sample_TypeSpermosphere:Time8                               -3.084
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control   7.229
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control  -5.865
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control    6.536
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith             1.439
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith             1.232
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith              1.893
    ## Sample_TypeSpermosphere:Time12:Location2Michigan            -0.839
    ## Sample_TypeSpermosphere:Time16:Location2Michigan            -0.613
    ## Sample_TypeSpermosphere:Time8:Location2Michigan              0.922
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota         2.123
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota         1.870
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota          1.961
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania         1.815
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania         0.995
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania          1.927
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                0.221
    ## Sample_TypeSpermosphere:Time16:Location2TVREC               -0.175
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                 1.472
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass            1.679
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass            0.888
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass             1.259
    ##                                                                        Pr(>|z|)
    ## (Intercept)                                                < 0.0000000000000002
    ## Sample_TypeSpermosphere                                                0.994631
    ## Time12                                                                 0.501874
    ## Time16                                                                 0.204319
    ## Time8                                                                  0.824556
    ## Location2Autoclaved_Control                                < 0.0000000000000002
    ## Location2EV_Smith                                                      0.026930
    ## Location2Michigan                                                      0.689843
    ## Location2North_Dakota                                                  0.000213
    ## Location2Pennsylvania                                                  0.440853
    ## Location2TVREC                                                         0.940481
    ## Location2Wiregrass                                                     0.684231
    ## Time12:Location2Autoclaved_Control                                     0.002281
    ## Time16:Location2Autoclaved_Control                            0.000003303125771
    ## Time8:Location2Autoclaved_Control                                      0.378588
    ## Time12:Location2EV_Smith                                               0.597922
    ## Time16:Location2EV_Smith                                               0.584657
    ## Time8:Location2EV_Smith                                                0.888682
    ## Time12:Location2Michigan                                               0.819638
    ## Time16:Location2Michigan                                               0.224364
    ## Time8:Location2Michigan                                                0.854263
    ## Time12:Location2North_Dakota                                           0.659453
    ## Time16:Location2North_Dakota                                           0.486948
    ## Time8:Location2North_Dakota                                            0.682514
    ## Time12:Location2Pennsylvania                                           0.513223
    ## Time16:Location2Pennsylvania                                           0.467036
    ## Time8:Location2Pennsylvania                                            0.839425
    ## Time12:Location2TVREC                                                  0.647547
    ## Time16:Location2TVREC                                                  0.288694
    ## Time8:Location2TVREC                                                   0.605072
    ## Time12:Location2Wiregrass                                              0.474099
    ## Time16:Location2Wiregrass                                              0.211842
    ## Time8:Location2Wiregrass                                               0.702051
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control                    0.006178
    ## Sample_TypeSpermosphere:Location2EV_Smith                              0.567441
    ## Sample_TypeSpermosphere:Location2Michigan                              0.749960
    ## Sample_TypeSpermosphere:Location2North_Dakota                          0.715506
    ## Sample_TypeSpermosphere:Location2Pennsylvania                          0.398874
    ## Sample_TypeSpermosphere:Location2TVREC                                 0.790674
    ## Sample_TypeSpermosphere:Location2Wiregrass                             0.415764
    ## Sample_TypeSpermosphere:Time12                                         0.000226
    ## Sample_TypeSpermosphere:Time16                                         0.003850
    ## Sample_TypeSpermosphere:Time8                                          0.002043
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control    0.000000000000487
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control    0.000000004479518
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control     0.000000000063074
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith                       0.150212
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith                       0.218072
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith                        0.058307
    ## Sample_TypeSpermosphere:Time12:Location2Michigan                       0.401557
    ## Sample_TypeSpermosphere:Time16:Location2Michigan                       0.540031
    ## Sample_TypeSpermosphere:Time8:Location2Michigan                        0.356540
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota                   0.033730
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota                   0.061531
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota                    0.049881
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania                   0.069546
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania                   0.319621
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania                    0.054014
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                          0.825382
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                          0.861016
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                           0.141073
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass                      0.093146
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass                      0.374438
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass                       0.208128
    ##                                                               
    ## (Intercept)                                                ***
    ## Sample_TypeSpermosphere                                       
    ## Time12                                                        
    ## Time16                                                        
    ## Time8                                                         
    ## Location2Autoclaved_Control                                ***
    ## Location2EV_Smith                                          *  
    ## Location2Michigan                                             
    ## Location2North_Dakota                                      ***
    ## Location2Pennsylvania                                         
    ## Location2TVREC                                                
    ## Location2Wiregrass                                            
    ## Time12:Location2Autoclaved_Control                         ** 
    ## Time16:Location2Autoclaved_Control                         ***
    ## Time8:Location2Autoclaved_Control                             
    ## Time12:Location2EV_Smith                                      
    ## Time16:Location2EV_Smith                                      
    ## Time8:Location2EV_Smith                                       
    ## Time12:Location2Michigan                                      
    ## Time16:Location2Michigan                                      
    ## Time8:Location2Michigan                                       
    ## Time12:Location2North_Dakota                                  
    ## Time16:Location2North_Dakota                                  
    ## Time8:Location2North_Dakota                                   
    ## Time12:Location2Pennsylvania                                  
    ## Time16:Location2Pennsylvania                                  
    ## Time8:Location2Pennsylvania                                   
    ## Time12:Location2TVREC                                         
    ## Time16:Location2TVREC                                         
    ## Time8:Location2TVREC                                          
    ## Time12:Location2Wiregrass                                     
    ## Time16:Location2Wiregrass                                     
    ## Time8:Location2Wiregrass                                      
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control        ** 
    ## Sample_TypeSpermosphere:Location2EV_Smith                     
    ## Sample_TypeSpermosphere:Location2Michigan                     
    ## Sample_TypeSpermosphere:Location2North_Dakota                 
    ## Sample_TypeSpermosphere:Location2Pennsylvania                 
    ## Sample_TypeSpermosphere:Location2TVREC                        
    ## Sample_TypeSpermosphere:Location2Wiregrass                    
    ## Sample_TypeSpermosphere:Time12                             ***
    ## Sample_TypeSpermosphere:Time16                             ** 
    ## Sample_TypeSpermosphere:Time8                              ** 
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control  ***
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith              
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith              
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith            .  
    ## Sample_TypeSpermosphere:Time12:Location2Michigan              
    ## Sample_TypeSpermosphere:Time16:Location2Michigan              
    ## Sample_TypeSpermosphere:Time8:Location2Michigan               
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota       *  
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota       .  
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota        *  
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania       .  
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania          
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania        .  
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                  
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass          .  
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass             
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass              
    ## 
    ## Phi coefficients (precision model with identity link):
    ##       Estimate Std. Error z value            Pr(>|z|)    
    ## (phi)   115.19       8.29   13.89 <0.0000000000000002 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Type of estimator: ML (maximum likelihood)
    ## Log-likelihood: 712.1 on 65 Df
    ## Pseudo R-squared: 0.9136
    ## Number of iterations: 75 (BFGS) + 2 (Fisher scoring)

    # Mean separation
    # Calculate estimated marginal means (least-squares means) for the treatment effect within combinations of type and time
    lsmeans <- emmeans(EvenStat.B, ~Location2 | Sample_Type * Time)

    # Perform multiple comparisons with Tukey adjustment
    Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE, Letters = letters)

    # Calculate observed evenness statistics
    ObservedEveness <- BAlphaCombined %>%
      group_by(Location2, Time, Sample_Type) %>%   
      summarise(MeanEvenness = mean(even), STDEvenness = sd(even)) %>% 
      arrange(-MeanEvenness)

    ## `summarise()` has grouped output by 'Location2', 'Time'. You can override using
    ## the `.groups` argument.

## Faith’s Phylogenetic Diversity Statistics

    # Fit a negative binomial regression model to the 'PD' (Faiths Phylogenetic Diversity) variable using the specified predictors
    RichStat.C <- glm.nb(PD ~ Sample_Type + Time + Location2 + Location2 * Time * Sample_Type, data = BAlphaCombined) 

    # Perform an ANOVA on the fitted model
    anova(RichStat.C)

    ## Analysis of Deviance Table
    ## 
    ## Model: Negative Binomial(268.5727), link: log
    ## 
    ## Response: PD
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ## 
    ##                            Df Deviance Resid. Df Resid. Dev
    ## NULL                                         383     5321.2
    ## Sample_Type                 1      1.9       382     5319.3
    ## Time                        3      5.1       379     5314.2
    ## Location2                   7   4174.8       372     1139.4
    ## Time:Location2             21     93.5       351     1045.9
    ## Sample_Type:Location2       7    383.1       344      662.8
    ## Sample_Type:Time            3      4.5       341      658.3
    ## Sample_Type:Time:Location2 21    146.9       320      511.3
    ##                                         Pr(>Chi)    
    ## NULL                                                
    ## Sample_Type                               0.1728    
    ## Time                                      0.1650    
    ## Location2                  < 0.00000000000000022 ***
    ## Time:Location2                  0.00000000003966 ***
    ## Sample_Type:Location2      < 0.00000000000000022 ***
    ## Sample_Type:Time                          0.2120    
    ## Sample_Type:Time:Location2 < 0.00000000000000022 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # Display a summary of the fitted model
    summary(RichStat.C)

    ## 
    ## Call:
    ## glm.nb(formula = PD ~ Sample_Type + Time + Location2 + Location2 * 
    ##     Time * Sample_Type, data = BAlphaCombined, init.theta = 268.5726887, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                                             Estimate Std. Error
    ## (Intercept)                                                 5.982650   0.032263
    ## Sample_TypeSpermosphere                                    -0.026493   0.045751
    ## Time12                                                      0.026690   0.045506
    ## Time16                                                      0.043367   0.045431
    ## Time8                                                      -0.056877   0.045896
    ## Location2Autoclaved_Control                                -0.963055   0.052561
    ## Location2EV_Smith                                           0.142046   0.045013
    ## Location2Michigan                                           0.095529   0.045206
    ## Location2North_Dakota                                      -0.274347   0.047059
    ## Location2Pennsylvania                                       0.062260   0.045348
    ## Location2TVREC                                              0.035202   0.045468
    ## Location2Wiregrass                                          0.120718   0.045100
    ## Time12:Location2Autoclaved_Control                         -0.110070   0.074900
    ## Time16:Location2Autoclaved_Control                         -0.314464   0.076488
    ## Time8:Location2Autoclaved_Control                          -0.176463   0.076416
    ## Time12:Location2EV_Smith                                    0.004421   0.063483
    ## Time16:Location2EV_Smith                                   -0.028607   0.063476
    ## Time8:Location2EV_Smith                                     0.079949   0.063786
    ## Time12:Location2Michigan                                   -0.035435   0.063870
    ## Time16:Location2Michigan                                   -0.061146   0.063844
    ## Time8:Location2Michigan                                     0.052079   0.064137
    ## Time12:Location2North_Dakota                               -0.052187   0.066575
    ## Time16:Location2North_Dakota                               -0.068470   0.066523
    ## Time8:Location2North_Dakota                                -0.027734   0.067101
    ## Time12:Location2Pennsylvania                               -0.053879   0.064131
    ## Time16:Location2Pennsylvania                               -0.052129   0.064020
    ## Time8:Location2Pennsylvania                                 0.031898   0.064401
    ## Time12:Location2TVREC                                       0.016468   0.064081
    ## Time16:Location2TVREC                                      -0.058948   0.064212
    ## Time8:Location2TVREC                                        0.115304   0.064313
    ## Time12:Location2Wiregrass                                  -0.010375   0.063647
    ## Time16:Location2Wiregrass                                  -0.019440   0.063572
    ## Time8:Location2Wiregrass                                    0.055158   0.063979
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control        -0.019815   0.074759
    ## Sample_TypeSpermosphere:Location2EV_Smith                   0.040223   0.063708
    ## Sample_TypeSpermosphere:Location2Michigan                   0.026487   0.064019
    ## Sample_TypeSpermosphere:Location2North_Dakota               0.004812   0.066727
    ## Sample_TypeSpermosphere:Location2Pennsylvania               0.037006   0.064188
    ## Sample_TypeSpermosphere:Location2TVREC                      0.036148   0.064358
    ## Sample_TypeSpermosphere:Location2Wiregrass                  0.049211   0.063804
    ## Sample_TypeSpermosphere:Time12                             -0.090876   0.064837
    ## Sample_TypeSpermosphere:Time16                             -0.118396   0.064823
    ## Sample_TypeSpermosphere:Time8                              -0.038871   0.065224
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control  0.854838   0.103545
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control  0.870077   0.105339
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   0.838271   0.105055
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith            0.063601   0.090124
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith            0.093354   0.090175
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith            -0.024702   0.090509
    ## Sample_TypeSpermosphere:Time12:Location2Michigan            0.011494   0.090845
    ## Sample_TypeSpermosphere:Time16:Location2Michigan            0.086950   0.090767
    ## Sample_TypeSpermosphere:Time8:Location2Michigan             0.020194   0.090970
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota        0.163021   0.094398
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota        0.102063   0.094650
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota         0.146423   0.094915
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania        0.091786   0.090989
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania        0.026200   0.091110
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania         0.061035   0.091210
    ## Sample_TypeSpermosphere:Time12:Location2TVREC               0.031205   0.091055
    ## Sample_TypeSpermosphere:Time16:Location2TVREC               0.113026   0.091185
    ## Sample_TypeSpermosphere:Time8:Location2TVREC               -0.086014   0.091414
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass           0.060983   0.090324
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass           0.070406   0.090321
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass            0.046396   0.090601
    ##                                                            z value
    ## (Intercept)                                                185.432
    ## Sample_TypeSpermosphere                                     -0.579
    ## Time12                                                       0.587
    ## Time16                                                       0.955
    ## Time8                                                       -1.239
    ## Location2Autoclaved_Control                                -18.323
    ## Location2EV_Smith                                            3.156
    ## Location2Michigan                                            2.113
    ## Location2North_Dakota                                       -5.830
    ## Location2Pennsylvania                                        1.373
    ## Location2TVREC                                               0.774
    ## Location2Wiregrass                                           2.677
    ## Time12:Location2Autoclaved_Control                          -1.470
    ## Time16:Location2Autoclaved_Control                          -4.111
    ## Time8:Location2Autoclaved_Control                           -2.309
    ## Time12:Location2EV_Smith                                     0.070
    ## Time16:Location2EV_Smith                                    -0.451
    ## Time8:Location2EV_Smith                                      1.253
    ## Time12:Location2Michigan                                    -0.555
    ## Time16:Location2Michigan                                    -0.958
    ## Time8:Location2Michigan                                      0.812
    ## Time12:Location2North_Dakota                                -0.784
    ## Time16:Location2North_Dakota                                -1.029
    ## Time8:Location2North_Dakota                                 -0.413
    ## Time12:Location2Pennsylvania                                -0.840
    ## Time16:Location2Pennsylvania                                -0.814
    ## Time8:Location2Pennsylvania                                  0.495
    ## Time12:Location2TVREC                                        0.257
    ## Time16:Location2TVREC                                       -0.918
    ## Time8:Location2TVREC                                         1.793
    ## Time12:Location2Wiregrass                                   -0.163
    ## Time16:Location2Wiregrass                                   -0.306
    ## Time8:Location2Wiregrass                                     0.862
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control         -0.265
    ## Sample_TypeSpermosphere:Location2EV_Smith                    0.631
    ## Sample_TypeSpermosphere:Location2Michigan                    0.414
    ## Sample_TypeSpermosphere:Location2North_Dakota                0.072
    ## Sample_TypeSpermosphere:Location2Pennsylvania                0.577
    ## Sample_TypeSpermosphere:Location2TVREC                       0.562
    ## Sample_TypeSpermosphere:Location2Wiregrass                   0.771
    ## Sample_TypeSpermosphere:Time12                              -1.402
    ## Sample_TypeSpermosphere:Time16                              -1.826
    ## Sample_TypeSpermosphere:Time8                               -0.596
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control   8.256
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control   8.260
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control    7.979
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith             0.706
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith             1.035
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith             -0.273
    ## Sample_TypeSpermosphere:Time12:Location2Michigan             0.127
    ## Sample_TypeSpermosphere:Time16:Location2Michigan             0.958
    ## Sample_TypeSpermosphere:Time8:Location2Michigan              0.222
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota         1.727
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota         1.078
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota          1.543
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania         1.009
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania         0.288
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania          0.669
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                0.343
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                1.240
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                -0.941
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass            0.675
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass            0.780
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass             0.512
    ##                                                                        Pr(>|z|)
    ## (Intercept)                                                < 0.0000000000000002
    ## Sample_TypeSpermosphere                                                 0.56254
    ## Time12                                                                  0.55753
    ## Time16                                                                  0.33980
    ## Time8                                                                   0.21525
    ## Location2Autoclaved_Control                                < 0.0000000000000002
    ## Location2EV_Smith                                                       0.00160
    ## Location2Michigan                                                       0.03458
    ## Location2North_Dakota                                       0.00000000554726116
    ## Location2Pennsylvania                                                   0.16977
    ## Location2TVREC                                                          0.43880
    ## Location2Wiregrass                                                      0.00744
    ## Time12:Location2Autoclaved_Control                                      0.14168
    ## Time16:Location2Autoclaved_Control                          0.00003934919421709
    ## Time8:Location2Autoclaved_Control                                       0.02093
    ## Time12:Location2EV_Smith                                                0.94448
    ## Time16:Location2EV_Smith                                                0.65222
    ## Time8:Location2EV_Smith                                                 0.21006
    ## Time12:Location2Michigan                                                0.57904
    ## Time16:Location2Michigan                                                0.33820
    ## Time8:Location2Michigan                                                 0.41679
    ## Time12:Location2North_Dakota                                            0.43311
    ## Time16:Location2North_Dakota                                            0.30335
    ## Time8:Location2North_Dakota                                             0.67937
    ## Time12:Location2Pennsylvania                                            0.40083
    ## Time16:Location2Pennsylvania                                            0.41549
    ## Time8:Location2Pennsylvania                                             0.62038
    ## Time12:Location2TVREC                                                   0.79719
    ## Time16:Location2TVREC                                                   0.35860
    ## Time8:Location2TVREC                                                    0.07300
    ## Time12:Location2Wiregrass                                               0.87051
    ## Time16:Location2Wiregrass                                               0.75977
    ## Time8:Location2Wiregrass                                                0.38862
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control                     0.79097
    ## Sample_TypeSpermosphere:Location2EV_Smith                               0.52780
    ## Sample_TypeSpermosphere:Location2Michigan                               0.67906
    ## Sample_TypeSpermosphere:Location2North_Dakota                           0.94251
    ## Sample_TypeSpermosphere:Location2Pennsylvania                           0.56426
    ## Sample_TypeSpermosphere:Location2TVREC                                  0.57434
    ## Sample_TypeSpermosphere:Location2Wiregrass                              0.44054
    ## Sample_TypeSpermosphere:Time12                                          0.16103
    ## Sample_TypeSpermosphere:Time16                                          0.06778
    ## Sample_TypeSpermosphere:Time8                                           0.55120
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control < 0.0000000000000002
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control   0.00000000000000147
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith                        0.48037
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith                        0.30055
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith                         0.78491
    ## Sample_TypeSpermosphere:Time12:Location2Michigan                        0.89932
    ## Sample_TypeSpermosphere:Time16:Location2Michigan                        0.33809
    ## Sample_TypeSpermosphere:Time8:Location2Michigan                         0.82432
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota                    0.08418
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota                    0.28089
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota                     0.12291
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania                    0.31309
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania                    0.77368
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania                     0.50339
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                           0.73182
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                           0.21515
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                            0.34674
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass                       0.49957
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass                       0.43568
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass                        0.60859
    ##                                                               
    ## (Intercept)                                                ***
    ## Sample_TypeSpermosphere                                       
    ## Time12                                                        
    ## Time16                                                        
    ## Time8                                                         
    ## Location2Autoclaved_Control                                ***
    ## Location2EV_Smith                                          ** 
    ## Location2Michigan                                          *  
    ## Location2North_Dakota                                      ***
    ## Location2Pennsylvania                                         
    ## Location2TVREC                                                
    ## Location2Wiregrass                                         ** 
    ## Time12:Location2Autoclaved_Control                            
    ## Time16:Location2Autoclaved_Control                         ***
    ## Time8:Location2Autoclaved_Control                          *  
    ## Time12:Location2EV_Smith                                      
    ## Time16:Location2EV_Smith                                      
    ## Time8:Location2EV_Smith                                       
    ## Time12:Location2Michigan                                      
    ## Time16:Location2Michigan                                      
    ## Time8:Location2Michigan                                       
    ## Time12:Location2North_Dakota                                  
    ## Time16:Location2North_Dakota                                  
    ## Time8:Location2North_Dakota                                   
    ## Time12:Location2Pennsylvania                                  
    ## Time16:Location2Pennsylvania                                  
    ## Time8:Location2Pennsylvania                                   
    ## Time12:Location2TVREC                                         
    ## Time16:Location2TVREC                                         
    ## Time8:Location2TVREC                                       .  
    ## Time12:Location2Wiregrass                                     
    ## Time16:Location2Wiregrass                                     
    ## Time8:Location2Wiregrass                                      
    ## Sample_TypeSpermosphere:Location2Autoclaved_Control           
    ## Sample_TypeSpermosphere:Location2EV_Smith                     
    ## Sample_TypeSpermosphere:Location2Michigan                     
    ## Sample_TypeSpermosphere:Location2North_Dakota                 
    ## Sample_TypeSpermosphere:Location2Pennsylvania                 
    ## Sample_TypeSpermosphere:Location2TVREC                        
    ## Sample_TypeSpermosphere:Location2Wiregrass                    
    ## Sample_TypeSpermosphere:Time12                                
    ## Sample_TypeSpermosphere:Time16                             .  
    ## Sample_TypeSpermosphere:Time8                                 
    ## Sample_TypeSpermosphere:Time12:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time16:Location2Autoclaved_Control ***
    ## Sample_TypeSpermosphere:Time8:Location2Autoclaved_Control  ***
    ## Sample_TypeSpermosphere:Time12:Location2EV_Smith              
    ## Sample_TypeSpermosphere:Time16:Location2EV_Smith              
    ## Sample_TypeSpermosphere:Time8:Location2EV_Smith               
    ## Sample_TypeSpermosphere:Time12:Location2Michigan              
    ## Sample_TypeSpermosphere:Time16:Location2Michigan              
    ## Sample_TypeSpermosphere:Time8:Location2Michigan               
    ## Sample_TypeSpermosphere:Time12:Location2North_Dakota       .  
    ## Sample_TypeSpermosphere:Time16:Location2North_Dakota          
    ## Sample_TypeSpermosphere:Time8:Location2North_Dakota           
    ## Sample_TypeSpermosphere:Time12:Location2Pennsylvania          
    ## Sample_TypeSpermosphere:Time16:Location2Pennsylvania          
    ## Sample_TypeSpermosphere:Time8:Location2Pennsylvania           
    ## Sample_TypeSpermosphere:Time12:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time16:Location2TVREC                 
    ## Sample_TypeSpermosphere:Time8:Location2TVREC                  
    ## Sample_TypeSpermosphere:Time12:Location2Wiregrass             
    ## Sample_TypeSpermosphere:Time16:Location2Wiregrass             
    ## Sample_TypeSpermosphere:Time8:Location2Wiregrass              
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(268.5727) family taken to be 1)
    ## 
    ##     Null deviance: 5321.16  on 383  degrees of freedom
    ## Residual deviance:  511.35  on 320  degrees of freedom
    ## AIC: 3938.8
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  268.6 
    ##           Std. Err.:  45.0 
    ## 
    ##  2 x log-likelihood:  -3808.77

    # Calculate the ratio of deviance to degrees of freedom for the fitted model
    RichStat.C$deviance / RichStat.C$df.residual

    ## [1] 1.597959

    # Mean separation
    # Calculate estimated marginal means (least-squares means) for the treatment effect within combinations of type and time
    lsmeans <- emmeans(RichStat.C, ~Location2 | Sample_Type * Time)

    # Perform multiple comparisons with Tukey adjustment and display results
    Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE, Letters = letters)

    # Calculate observed Faith's Phylogenetic Diversity statistics
    ObservedFPD <- BAlphaCombined %>%
      group_by(Location2, Time, Sample_Type) %>%   
      summarise(MeanFPD = mean(PD), STFPD = sd(PD)) %>% 
      arrange(-MeanFPD)

    ## `summarise()` has grouped output by 'Location2', 'Time'. You can override using
    ## the `.groups` argument.

\####Conjoined Figure

    figure1b <- ggpubr::ggarrange(bac.richness,
                                           bac.even,
                                           bac.phylo.div,
                                           labels = "auto",
                                           nrow = 3, ncol = 1, common.legend = T)

    ## Warning: Removed 648 rows containing missing values or values outside the scale range
    ## (`geom_pwc()`).
    ## Removed 648 rows containing missing values or values outside the scale range
    ## (`geom_pwc()`).

    #Arranging the figures
    #Note since figures are so large this is unnecessary

# Beta Diversity

## Beta-Diversity using Weighted Unifrac Only Bacteria

    ##### All Samples ########
    bac.unifrac = UniFrac(BSpermosphere_RDS_CSS, weighted = TRUE) # create weighted unifrac distance matrix

    ## Warning in UniFrac(BSpermosphere_RDS_CSS, weighted = TRUE): Randomly assigning
    ## root as -- BOTU_7931 -- in the phylogenetic tree in the data you provided.

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)  # Set seed for reproducibility
    unifracstats <- adonis2(bac.unifrac ~ Location_EC * Sample_Type * Time, as(sample_data(BSpermosphere_RDS_CSS), "data.frame"))  # Perform PERMANOVA
    #the model is significant therefore we will subset by time.

    # Perform ordination using MDS with weighted UniFrac distance
    global.ord.uni <- phyloseq::ordinate(BSpermosphere_RDS_CSS, "MDS", "unifrac", weighted = TRUE)

    ## Warning in UniFrac(physeq, ...): Randomly assigning root as -- BOTU_5859 -- in
    ## the phylogenetic tree in the data you provided.

    # Calculate percentage variation on axis one
    variation.axis1 <- round(100 * global.ord.uni$values$Relative_eig[[1]], 2)

    # Calculate percentage variation on axis two
    variation.axis2 <- round(100 * global.ord.uni$values$Relative_eig[[2]], 2)

    # Extract ordination vectors
    bac.ggplot.data.uni <- data.frame(global.ord.uni$vectors)
    bac.ggplot.data.uni$Sample_ID <- rownames(bac.ggplot.data.uni)

    # Merge ordination vectors with sample metadata
    bac.ggplot.data.uni2 <- left_join(data.frame(BSpermosphere_RDS_CSS@sam_data), bac.ggplot.data.uni, by = "Sample_ID")

    # Order the time points in chronological order
    bac.ggplot.data.uni2$Time.Point <- factor(bac.ggplot.data.uni2$Time, levels = c("0", "8", "12", "16"))

    # Create a ggplot for weighted UniFrac distance

    global.uni <- ggplot(bac.ggplot.data.uni2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type, color = Location_EC)) +
      geom_point(size = 4, alpha = 0.7) +  # Add points with specified size and transparency
      scale_shape_manual(values = c(21, 22, 24, 25)) +  # Customize shapes for different types
      scale_fill_manual(values = cbbPalette) +  # Customize fill colors for different treatments
      scale_color_manual(values = cbbPalette) +  # Customize colors for different treatments
      xlab(paste("PcoA1 -", variation.axis1, "%")) +  # Label for x-axis with percentage variation 1
      ylab(paste("PcoA2 -", variation.axis2, "%")) +  # Label for y-axis with percentage variation 2
      theme(legend.position = "bottom", plot.title = element_text(size = 10, face = "bold")) +  # Customize theme
      ggtitle("Weighted Unifrac")  # Add title to the plot
    global.uni  # Display the plot

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Global%20Weighted%20Unifrac-1.png)
\### 0 hours; Unifrac \#####

    # Filter to time point
    bac_sperm_0 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "0") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)  # Set seed for reproducibility
    prok.dist.uni.0hrs = UniFrac(bac_sperm_0, weighted = TRUE)  # Create distance matrix

    ## Warning in UniFrac(bac_sperm_0, weighted = TRUE): Randomly assigning root as --
    ## BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.0.uni <- adonis2(prok.dist.uni.0hrs ~ Location_EC * Sample_Type, as(sample_data(bac_sperm_0), "data.frame"))  # Perform PERMANOVA to test for significant changes
    adonis.0.uni  # Display PERMANOVA results

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.0hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_0), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15  2.47402 0.90938 53.523  0.001 ***
    ## Residual 80  0.24653 0.09062                  
    ## Total    95  2.72054 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the p-value for the factor
    pvalue <- adonis.0.uni$`Pr(>F)`[[1]]
    R2 <- adonis.0.uni$R2[[1]]

    # Plotting
    ordination.unifrac.pcoa.0hrs <- ordinate(bac_sperm_0, "PCoA", distance = prok.dist.uni.0hrs)  # Calculate the resemblance and ordinate using PCoA
    uni.0.values <- data.frame(ordination.unifrac.pcoa.0hrs$values)

    # Percent variation
    variation.axis1.0.uni <- round(100 * uni.0.values$Relative_eig[[1]], 2)
    variation.axis2.0.uni <- round(100 * uni.0.values$Relative_eig[[2]], 2)

    # Data from ordinations
    uni.0.vectors <- data.frame(ordination.unifrac.pcoa.0hrs$vectors)  # Positions of your points on the PCoA graph
    uni.0.vectors$Sample_ID <- rownames(uni.0.vectors)
    uni.0.vectors2 <- left_join(data.frame(bac_sperm_0@sam_data), uni.0.vectors, by = "Sample_ID")

    # Create a ggplot for UniFrac distance at 0 hours
    pcoa.unifrac.0hrs.new <- ggplot(uni.0.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size = 4, alpha = 0.7) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.0.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.0.uni, "%")) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) +
      ggtitle(paste("Planting")) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
    pcoa.unifrac.0hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/0hours%20Weighted%20Unifrac-1.png)

    # Dispersion for each Time-Point
    beta.disp.0.uni <- betadisper(prok.dist.uni.0hrs, bac_sperm_0@sam_data$Location_EC)
    permutest(beta.disp.0.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 0.10714 0.0153060 47.781    999  0.001 ***
    ## Residuals 88 0.02819 0.0003203                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.0.unib <- betadisper(prok.dist.uni.0hrs, bac_sperm_0@sam_data$Sample_Type)
    permutest(beta.disp.0.unib)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00262 0.0026151 0.3202    999  0.593
    ## Residuals 94 0.76770 0.0081671

### 8 hours; Unifrac

    # Filter to time point
    bac_sperm_8 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "8") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)  # Set seed for reproducibility
    prok.dist.uni.8hrs = UniFrac(bac_sperm_8, weighted = TRUE)  # Create distance matrix

    ## Warning in UniFrac(bac_sperm_8, weighted = TRUE): Randomly assigning root as --
    ## BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.8.uni <- adonis2(prok.dist.uni.8hrs ~ Location_EC * Sample_Type, as(sample_data(bac_sperm_8), "data.frame"))  # Perform PERMANOVA to test for significant changes
    adonis.8.uni  # Display PERMANOVA results

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.8hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_8), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15   2.2455 0.88114 39.538  0.001 ***
    ## Residual 80   0.3029 0.11886                  
    ## Total    95   2.5484 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the p-value for the factor
    pvalue <- adonis.8.uni$`Pr(>F)`[[1]]
    R2 <- adonis.8.uni$R2[[1]]

    # Plotting
    ordination.unifrac.pcoa.8hrs <- ordinate(bac_sperm_8, "PCoA", distance = prok.dist.uni.8hrs)  # Calculate the resemblance and ordinate using PCoA
    uni.8.values <- data.frame(ordination.unifrac.pcoa.8hrs$values)

    # Percent variation
    variation.axis1.8.uni <- round(100 * uni.8.values$Relative_eig[[1]], 2)
    variation.axis2.8.uni <- round(100 * uni.8.values$Relative_eig[[2]], 2)

    # Data from ordinations
    uni.8.vectors <- data.frame(ordination.unifrac.pcoa.8hrs$vectors)  # Positions of your points on the PCoA graph
    uni.8.vectors$Sample_ID <- rownames(uni.8.vectors)
    uni.8.vectors2 <- left_join(data.frame(bac_sperm_8@sam_data), uni.8.vectors, by = "Sample_ID")

    # Create a ggplot for UniFrac distance at 8 hours
    pcoa.unifrac.8hrs.new <- ggplot(uni.8.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size = 4, alpha = 0.7) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.8.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.8.uni, "%")) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) +
      ggtitle(paste("8 Hours Postplanting")) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
    pcoa.unifrac.8hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/8hours%20Weighted%20Unifrac-1.png)

    # Dispersion for each Time-Point
    beta.disp.8.uni <- betadisper(prok.dist.uni.8hrs, bac_sperm_8@sam_data$Location_EC)
    permutest(beta.disp.8.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 0.112057 0.0160082 18.808    999  0.001 ***
    ## Residuals 88 0.074901 0.0008511                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.8.unib <- betadisper(prok.dist.uni.8hrs, bac_sperm_8@sam_data$Sample_Type)
    permutest(beta.disp.8.unib)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.01060 0.0106041 1.7916    999  0.195
    ## Residuals 94 0.55637 0.0059188

### 12; Unifrac

    # Filter to time point
    bac_sperm_12 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "12") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)  # Set seed for reproducibility
    prok.dist.uni.12hrs = UniFrac(bac_sperm_12, weighted = TRUE)  # Create distance matrix

    ## Warning in UniFrac(bac_sperm_12, weighted = TRUE): Randomly assigning root as
    ## -- BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.12.uni <- adonis2(prok.dist.uni.12hrs ~ Location_EC * Sample_Type, as(sample_data(bac_sperm_12), "data.frame"))  # Perform PERMANOVA to test for significant changes
    adonis.12.uni  # Display PERMANOVA results

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.12hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_12), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15  2.25765 0.86355 33.752  0.001 ***
    ## Residual 80  0.35674 0.13645                  
    ## Total    95  2.61438 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the p-value for the factor
    pvalue <- adonis.12.uni$`Pr(>F)`[[1]]
    R2 <- adonis.12.uni$R2[[1]]

    # Plotting
    ordination.unifrac.pcoa.12hrs <- ordinate(bac_sperm_12, "PCoA", distance = prok.dist.uni.12hrs)  # Calculate the resemblance and ordinate using PCoA
    uni.12.values <- data.frame(ordination.unifrac.pcoa.12hrs$values)

    # Percent variation
    variation.axis1.12.uni <- round(100 * uni.12.values$Relative_eig[[1]], 2)
    variation.axis2.12.uni <- round(100 * uni.12.values$Relative_eig[[2]], 2)

    # Data from ordinations
    uni.12.vectors <- data.frame(ordination.unifrac.pcoa.12hrs$vectors)  # Positions of your points on the PCoA graph
    uni.12.vectors$Sample_ID <- rownames(uni.12.vectors)
    uni.12.vectors2 <- left_join(data.frame(bac_sperm_12@sam_data), uni.12.vectors, by = "Sample_ID")

    # Create a ggplot for UniFrac distance at 12 growth stage
    pcoa.unifrac.12hrs.new <- ggplot(uni.12.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size = 4, alpha = 0.7) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.12.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.12.uni, "%")) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) +
      ggtitle(paste("12 Growth Stage")) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
    pcoa.unifrac.12hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/12hours%20Weighted%20Unifrac-1.png)

    # Dispersion for each Time-Point
    beta.disp.12.uni <- betadisper(prok.dist.uni.12hrs, bac_sperm_12@sam_data$Location_EC)
    permutest(beta.disp.12.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 0.11065 0.0158074 8.6256    999  0.001 ***
    ## Residuals 88 0.16127 0.0018326                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.12.unib <- betadisper(prok.dist.uni.12hrs, bac_sperm_12@sam_data$Sample_Type)
    permutest(beta.disp.12.unib)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00400 0.0040028 0.8443    999  0.403
    ## Residuals 94 0.44565 0.0047409

### 16; Unifrac

    # Filter to time point
    bac_sperm_16 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "16") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)  # Set seed for reproducibility
    prok.dist.uni.16hrs = UniFrac(bac_sperm_16, weighted = TRUE)  # Create distance matrix

    ## Warning in UniFrac(bac_sperm_16, weighted = TRUE): Randomly assigning root as
    ## -- BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.16.uni <- adonis2(prok.dist.uni.16hrs ~ Location_EC * Sample_Type, as(sample_data(bac_sperm_16), "data.frame"))  # Perform PERMANOVA to test for significant changes
    adonis.16.uni  # Display PERMANOVA results

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.16hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_16), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15   2.5531 0.79735 20.984  0.001 ***
    ## Residual 80   0.6489 0.20265                  
    ## Total    95   3.2020 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the p-value for the factor
    pvalue <- adonis.16.uni$`Pr(>F)`[[1]]
    R2 <- adonis.16.uni$R2[[1]]

    # Plotting
    ordination.unifrac.pcoa.16hrs <- ordinate(bac_sperm_16, "PCoA", distance = prok.dist.uni.16hrs)  # Calculate the resemblance and ordinate using PCoA
    uni.16.values <- data.frame(ordination.unifrac.pcoa.16hrs$values)

    # Percent variation
    variation.axis1.16.uni <- round(100 * uni.16.values$Relative_eig[[1]], 2)
    variation.axis2.16.uni <- round(100 * uni.16.values$Relative_eig[[2]], 2)

    # Data from ordinations
    uni.16.vectors <- data.frame(ordination.unifrac.pcoa.16hrs$vectors)  # Positions of your points on the PCoA graph
    uni.16.vectors$Sample_ID <- rownames(uni.16.vectors)
    uni.16.vectors2 <- left_join(data.frame(bac_sperm_16@sam_data), uni.16.vectors, by = "Sample_ID")

    # Create a ggplot for UniFrac distance at 16 growth stage
    pcoa.unifrac.16hrs.new <- ggplot(uni.16.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size = 4, alpha = 0.7) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.16.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.16.uni, "%")) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) +
      ggtitle(paste("16 Growth Stage")) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
    pcoa.unifrac.16hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/16hours%20Weighted%20Unifrac-1.png)

    # Dispersion for each Time-Point
    beta.disp.16.uni <- betadisper(prok.dist.uni.16hrs, bac_sperm_16@sam_data$Location_EC)
    permutest(beta.disp.16.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 0.21188 0.0302690 10.333    999  0.001 ***
    ## Residuals 88 0.25778 0.0029294                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.16.unib <- betadisper(prok.dist.uni.16hrs, bac_sperm_16@sam_data$Sample_Type)
    permutest(beta.disp.16.unib)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00695 0.0069469 0.9885    999  0.336
    ## Residuals 94 0.66060 0.0070277

## Beta-Diversity using UnWeighted Unifrac Only Bacteria

    ##### All Samples ########
    # Create unweighted UniFrac distance matrix
    bac.unifrac = UniFrac(BSpermosphere_RDS_CSS, weighted = FALSE)

    ## Warning in UniFrac(BSpermosphere_RDS_CSS, weighted = FALSE): Randomly assigning
    ## root as -- BOTU_1548 -- in the phylogenetic tree in the data you provided.

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)
    adonis2(bac.unifrac ~ Location_EC * Sample_Type * Time, as(sample_data(BSpermosphere_RDS_CSS), "data.frame"))

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = bac.unifrac ~ Location_EC * Sample_Type * Time, data = as(sample_data(BSpermosphere_RDS_CSS), "data.frame"))
    ##           Df SumOfSqs     R2      F Pr(>F)    
    ## Model     63   19.066 0.5939 7.4283  0.001 ***
    ## Residual 320   13.037 0.4061                  
    ## Total    383   32.103 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # Perform ordination using MDS with unweighted UniFrac distance
    global.ord.uni <- phyloseq::ordinate(BSpermosphere_RDS_CSS, "MDS", "unifrac", weighted = FALSE)

    ## Warning in UniFrac(physeq, ...): Randomly assigning root as -- BOTU_5859 -- in
    ## the phylogenetic tree in the data you provided.

    # Calculate percentage variation on axis one
    variation.axis1 <- round(100 * global.ord.uni$values$Relative_eig[[1]], 2)

    # Calculate percentage variation on axis two
    variation.axis2 <- round(100 * global.ord.uni$values$Relative_eig[[2]], 2)

    # Extract ordination vectors
    bac.ggplot.data.uni <- data.frame(global.ord.uni$vectors)
    bac.ggplot.data.uni$Sample_ID <- rownames(bac.ggplot.data.uni)

    # Merge ordination vectors with sample metadata
    bac.ggplot.data.uni2 <- left_join(data.frame(BSpermosphere_RDS_CSS@sam_data), bac.ggplot.data.uni, by = "Sample_ID")

    # Order the time points in chronological order
    bac.ggplot.data.uni2$Time.Point <- factor(bac.ggplot.data.uni2$Time, levels = c("0", "8", "12","16"))

    # Create a ggplot for unweighted UniFrac distance
    global.uni.unweighted <- ggplot(bac.ggplot.data.uni2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type, color = Location_EC)) +
      geom_point(size = 4, alpha = 0.7) +  # Add points with specified size and transparency
      scale_shape_manual(values = c(21, 22, 24, 25)) +  # Customize shapes for different types
      scale_fill_manual(values = cbbPalette) +  # Customize fill colors for different treatments
      scale_color_manual(values = cbbPalette) +  # Customize colors for different treatments
      xlab(paste("PcoA1 -", variation.axis1, "%")) +  # Label for x-axis with percentage variation 1
      ylab(paste("PcoA2 -", variation.axis2, "%")) +  # Label for y-axis with percentage variation 2
      theme(legend.position = "bottom", plot.title = element_text(size = 10, face = "bold")) +  # Customize theme
      ggtitle("Unweighted Unifrac")  # Add title to the plot
    global.uni.unweighted  # Display the plot

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/Global%20Unweighted%20Unifrac-1.png)
\### 0 hours; Unweighted Unifrac \#####

    # Filter to time point
    bac_sperm_0 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "0") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)
    prok.dist.uni.0hrs = UniFrac(bac_sperm_0, weighted = FALSE)  # Create distance matrix

    ## Warning in UniFrac(bac_sperm_0, weighted = FALSE): Randomly assigning root as
    ## -- BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.0.uni <- adonis2(prok.dist.uni.0hrs ~ Location_EC * Sample_Type, as(sample_data(bac_sperm_0), "data.frame"))  # Perform PERMANOVA to test for significant changes
    adonis.0.uni  # Display PERMANOVA results

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.0hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_0), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15   5.0873 0.60766 8.2604  0.001 ***
    ## Residual 80   3.2846 0.39234                  
    ## Total    95   8.3718 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the p-value for the factor
    pvalue <- adonis.0.uni$`Pr(>F)`[[1]]
    R2 <- adonis.0.uni$R2[[1]]

    # Plotting
    ordination.unifrac.pcoa.0hrs <- ordinate(bac_sperm_0, "PCoA", distance = prok.dist.uni.0hrs)  # Calculate the resemblance and ordinate using PCoA
    uni.0.values <- data.frame(ordination.unifrac.pcoa.0hrs$values)

    # Percent variation
    variation.axis1.0.uni <- round(100 * uni.0.values$Relative_eig[[1]], 2)
    variation.axis2.0.uni <- round(100 * uni.0.values$Relative_eig[[2]], 2)

    # Data from ordinations
    uni.0.vectors <- data.frame(ordination.unifrac.pcoa.0hrs$vectors)  # Positions of your points on the PCoA graph
    uni.0.vectors$Sample_ID <- rownames(uni.0.vectors)
    uni.0.vectors2 <- left_join(data.frame(bac_sperm_0@sam_data), uni.0.vectors, by = "Sample_ID")

    # Create a ggplot for unweighted UniFrac distance at 0 hours
    un.pcoa.unifrac.0hrs.new <- ggplot(uni.0.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size = 4, alpha = 0.7) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.0.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.0.uni, "%")) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) +
      ggtitle(paste("Planting")) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.key = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))
    un.pcoa.unifrac.0hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/0hours%20Unweighted%20Unifrac-1.png)

    # Dispersion for each Time-Point
    beta.disp.0.uni <- betadisper(prok.dist.uni.0hrs, bac_sperm_0@sam_data$Location_EC)
    permutest(beta.disp.0.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 1.01238 0.14463 352.86    999  0.001 ***
    ## Residuals 88 0.03607 0.00041                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.0.uni.Ty <- betadisper(prok.dist.uni.0hrs, bac_sperm_0@sam_data$Sample_Type)

### 8 hours; Unweighted Unifrac

    bac_sperm_8 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "8") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)
    prok.dist.uni.8hrs = UniFrac(bac_sperm_8, weighted = F) # create distance matrix

    ## Warning in UniFrac(bac_sperm_8, weighted = F): Randomly assigning root as --
    ## BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.8.uni <- adonis2(prok.dist.uni.8hrs~Location_EC*Sample_Type, as(sample_data(bac_sperm_8), "data.frame")) #Are there significant changes?
    adonis.8.uni

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.8hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_8), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15   4.7877 0.59083 7.7012  0.001 ***
    ## Residual 80   3.3156 0.40917                  
    ## Total    95   8.1033 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the pvalue for the factor

    pvalue <- adonis.8.uni$`Pr(>F)`[[1]]
    R2 <- adonis.8.uni$R2[[1]]

    # plotting
    ordination.unifrac.pcoa.8hrs <- ordinate(bac_sperm_8, "PCoA", distance = prok.dist.uni.8hrs) # calculate the resemblance and ordinate using PCoA
    uni.8.values <- data.frame(ordination.unifrac.pcoa.8hrs$values)

    # percent variation
    variation.axis1.8.uni <- round(100*uni.8.values$Relative_eig[[1]], 2)
    variation.axis2.8.uni <- round(100*uni.8.values$Relative_eig[[2]], 2)

    # data from ordinations
    uni.8.vectors <- data.frame(ordination.unifrac.pcoa.8hrs$vectors) # positions of your points on the pcoa graph
    uni.8.vectors$Sample_ID <- rownames(uni.8.vectors)
    uni.8.vectors2 <- left_join(data.frame(bac_sperm_8@sam_data), uni.8.vectors,  by = "Sample_ID")
      
      un.pcoa.unifrac.8hrs.new <- ggplot(uni.8.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size=4, alpha = 0.7) +
      scale_shape_manual(values = c(21,24)) +
      scale_fill_manual(values=cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.8.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.8.uni, "%")) +
      guides(fill=guide_legend(override.aes=list(shape=21))) +
      ggtitle(paste("8 Hours Postplanting")) +
      theme_bw() +
      theme(legend.position="bottom",
            legend.title=element_blank(),
            legend.key = element_blank(),
            plot.title=element_text(size=10, face = "bold"))
    un.pcoa.unifrac.8hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/8hours%20Unweighted%20Unifrac-1.png)

    #Dispersion for each Time-Point
    beta.disp.8.uni <- betadisper(prok.dist.uni.8hrs, bac_sperm_8@sam_data$Location_EC)
    permutest(beta.disp.8.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 1.00410 0.143443 126.77    999  0.001 ***
    ## Residuals 88 0.09958 0.001132                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.8.uni.b <- betadisper(prok.dist.uni.8hrs, bac_sperm_8@sam_data$Sample_Type)
    permutest(beta.disp.8.uni.b)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.01684 0.016841 0.7073    999  0.425
    ## Residuals 94 2.23817 0.023810

### 12; Unweighted Unifrac

    bac_sperm_12 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "12") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)
    prok.dist.uni.12hrs = UniFrac(bac_sperm_12, weighted = F) # create distance matrix

    ## Warning in UniFrac(bac_sperm_12, weighted = F): Randomly assigning root as --
    ## BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.12.uni <- adonis2(prok.dist.uni.12hrs~Location_EC*Sample_Type, as(sample_data(bac_sperm_12), "data.frame")) #Are there significant changes?
    adonis.12.uni

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.12hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_12), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15   4.1899 0.58298 7.4559  0.001 ***
    ## Residual 80   2.9971 0.41702                  
    ## Total    95   7.1870 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the pvalue for the factor

    pvalue <- adonis.12.uni$`Pr(>F)`[[1]]
    R2 <- adonis.12.uni$R2[[1]]

    # ploting
    ordination.unifrac.pcoa.12hrs <- ordinate(bac_sperm_12, "PCoA", distance = prok.dist.uni.12hrs) # calculate the resemblance and ordinate using PCoA
    uni.12.values <- data.frame(ordination.unifrac.pcoa.12hrs$values)

    # percent variation
    variation.axis1.12.uni <- round(100*uni.12.values$Relative_eig[[1]], 2)
    variation.axis2.12.uni <- round(100*uni.12.values$Relative_eig[[2]], 2)

    # data from ordinations
    uni.12.vectors <- data.frame(ordination.unifrac.pcoa.12hrs$vectors) # positions of your points on the pcoa graph
    uni.12.vectors$Sample_ID <- rownames(uni.12.vectors)
    uni.12.vectors2 <- left_join(data.frame(bac_sperm_12@sam_data), uni.12.vectors,  by = "Sample_ID")

      un.pcoa.unifrac.12hrs.new <- ggplot(uni.12.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size=4, alpha = 0.7) +
      scale_shape_manual(values = c(21,24)) +
      scale_fill_manual(values=cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.12.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.12.uni, "%")) +
      guides(fill=guide_legend(override.aes=list(shape=21))) +
      ggtitle(paste("12 Hours Postplanting")) +
      theme_bw() +
      theme(legend.position="bottom",
            legend.title=element_blank(),
            legend.key = element_blank(),
            plot.title=element_text(size=10, face = "bold"))
    un.pcoa.unifrac.12hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/12hours%20Unweighted%20Unifrac-1.png)

    #Dispersion for each Time-Point
    beta.disp.12.uni <- betadisper(prok.dist.uni.12hrs, bac_sperm_12@sam_data$Location_EC)
    permutest(beta.disp.12.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 0.85299 0.121855 80.345    999  0.001 ***
    ## Residuals 88 0.13347 0.001517                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.12.uni.b <- betadisper(prok.dist.uni.12hrs, bac_sperm_12@sam_data$Sample_Type)
    permutest(beta.disp.12.uni.b)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
    ## Groups     1 0.01066 0.010664 0.547    999   0.48
    ## Residuals 94 1.83272 0.019497

### 16; Unweighted Unifrac

    bac_sperm_16 <- BSpermosphere_RDS_CSS %>% 
      subset_samples(Time == "16") %>% 
      phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

    # PERMANOVA - testing for differences in centroids
    set.seed(12345)
    prok.dist.uni.16hrs = UniFrac(bac_sperm_16, weighted = F) # create distance matrix

    ## Warning in UniFrac(bac_sperm_16, weighted = F): Randomly assigning root as --
    ## BOTU_109446 -- in the phylogenetic tree in the data you provided.

    adonis.16.uni <- adonis2(prok.dist.uni.16hrs~Location_EC*Sample_Type, as(sample_data(bac_sperm_16), "data.frame")) #Are there significant changes?
    adonis.16.uni

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = prok.dist.uni.16hrs ~ Location_EC * Sample_Type, data = as(sample_data(bac_sperm_16), "data.frame"))
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    15   4.8477 0.58529 7.5271  0.001 ***
    ## Residual 80   3.4348 0.41471                  
    ## Total    95   8.2825 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # The R2 is the variation due to different factors 
    # pvalue is the pvalue for the factor

    pvalue <- adonis.16.uni$`Pr(>F)`[[1]]
    R2 <- adonis.16.uni$R2[[1]]

    # ploting
    ordination.unifrac.pcoa.16hrs <- ordinate(bac_sperm_16, "PCoA", distance = prok.dist.uni.16hrs) # calculate the resemblance and ordinate using PCoA
    uni.16.values <- data.frame(ordination.unifrac.pcoa.16hrs$values)

    # percent variation
    variation.axis1.16.uni <- round(100*uni.16.values$Relative_eig[[1]], 2)
    variation.axis2.16.uni <- round(100*uni.16.values$Relative_eig[[2]], 2)

    # data from ordinations
    uni.16.vectors <- data.frame(ordination.unifrac.pcoa.16hrs$vectors) # positions of your points on the pcoa graph
    uni.16.vectors$Sample_ID <- rownames(uni.16.vectors)
    uni.16.vectors2 <- left_join(data.frame(bac_sperm_16@sam_data), uni.16.vectors,  by = "Sample_ID")

      un.pcoa.unifrac.16hrs.new <- ggplot(uni.16.vectors2, aes(x = Axis.1, y = Axis.2, fill = Location_EC, shape = Sample_Type)) +
      geom_point(size=4, alpha = 0.7) +
      scale_shape_manual(values = c(21,24)) +
      scale_fill_manual(values=cbbPalette) +
      xlab(paste("PcoA1 -", variation.axis1.16.uni, "%")) + 
      ylab(paste("PcoA2 -", variation.axis2.16.uni, "%")) +
      guides(fill=guide_legend(override.aes=list(shape=21))) +
      ggtitle(paste("16 Hours Postplanting")) +
      theme_bw() +
      theme(legend.position="bottom",
            legend.title=element_blank(),
            legend.key = element_blank(),
            plot.title=element_text(size=10, face = "bold"))
    un.pcoa.unifrac.16hrs.new

![](LL_Prokaryote_Diversity_Analysis_files/figure-markdown_strict/16hours%20Unweighted%20Unifrac-1.png)

    #Dispersion for each Time-Point
    beta.disp.16.uni <- betadisper(prok.dist.uni.16hrs, bac_sperm_16@sam_data$Location_EC)
    permutest(beta.disp.16.uni)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
    ## Groups     7 0.99498 0.142141 151.47    999  0.001 ***
    ## Residuals 88 0.08258 0.000938                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    beta.disp.16.uni.b <- betadisper(prok.dist.uni.16hrs, bac_sperm_16@sam_data$Sample_Type)
    permutest(beta.disp.16.uni.b)

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00201 0.0020055 0.0825    999  0.787
    ## Residuals 94 2.28523 0.0243110

## combining beta diversity figures

    UnifracCombined <- ggpubr::ggarrange(pcoa.unifrac.0hrs.new,
                                   pcoa.unifrac.8hrs.new,
                                   pcoa.unifrac.12hrs.new,
                                   pcoa.unifrac.16hrs.new,
                                   un.pcoa.unifrac.0hrs.new,
                                   un.pcoa.unifrac.8hrs.new,
                                   un.pcoa.unifrac.12hrs.new,
                                   un.pcoa.unifrac.16hrs.new,
                                   labels = "auto",
                                   nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")
    #A-D is Weighted Unifrac, E-H is unweighted Unifrac
    #Arranging and Combining the Beta Diversity Figures figures
