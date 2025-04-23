-   [Bacterial Data Preprocessing](#bacterial-data-preprocessing)
-   [Read in data](#read-in-data)
-   [Decontaminate](#decontaminate)
-   [Taxonomy filtering](#taxonomy-filtering)
-   [Mock Community analysis](#mock-community-analysis)
-   [Data filtering](#data-filtering)
-   [Read Depth Analysis](#read-depth-analysis)
-   [Rarefaction Analysis](#rarefaction-analysis)
-   [Metagenome CSS normalization](#metagenome-css-normalization)

# Bacterial Data Preprocessing

    ##Loading Libraries ##

    library(phyloseq)
    library(decontam)
    library(vegan)
    library(tidyverse)
    library(metagenomeSeq)
    library(ggplot2)
    library(ggpubr)
    library(Biostrings)
    library(microbiome)

    ##### Set global options #####
    # no scientific notation
    options(scipen=10000) 
    #Since RMD files moved the working directories needed to be specified

    # color blind pallets used throughout 
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
    tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

# Read in data

    # Metadata #
    samp_dat_bac <- read.csv("HPC_Output/Spermosphere_16S_Metadata.csv", na.strings = "na")

    rownames(samp_dat_bac) <- samp_dat_bac$Sample_ID #row names must match OTU table headers
    SAMP.bac <- phyloseq::sample_data(samp_dat_bac)

    #OTU table # remove hash and space in OTU_ID header manually
    #OTU TABLE is too large therefore it is included in a zip file, Github file size exceeded
    #otu_bac <- read.csv("HPC_Output/FASTA _and_OTU_TABLE/otu_table_16S_Bacteria.txt", sep = "\t") #txt file is tab separated
    #colnames(otu_bac) <- sub("_.*", "", colnames(otu_bac)) #removing _S*_merged
    #rownames(otu_bac) <- otu_bac$OTU
    #otu_bac <- otu_bac[,-1]
    #OTU.bac <- phyloseq::otu_table(otu_bac, taxa_are_rows = TRUE)
    #any(is.na(otu_bac)) # no NA in the OTU table
    ###IN ORDER TO KNIT THIS THE LINES ABOVE ARE COMMENTED OUT FILE IS COMPRESSED ###
    # Taxonomy #
    taxonomy.bac <- read.csv("HPC_Output/16s_taxonomy_SINTAX.txt", sep = "")
    rownames(taxonomy.bac) <- taxonomy.bac$OTU
    TAX.bac <- phyloseq::tax_table(as.matrix(taxonomy.bac))

    # Fasta # fasta is too large so not running for now
    # FASTA is too large therefore it is included in a zip file. Github file size exceeded.
    #FASTA.bac <- readDNAStringSet("HPC_Output/FASTA _and_OTU_TABLE/otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
    ###IN ORDER TO KNIT THIS THE LINE ABOVE IS COMMENTED OUT FILE IS COMPRESSED ###
    # Phylogentic tree #
    tree <- phyloseq::read_tree("HPC_Output/otu_tree.tre")

    ###### Create Initial Phyloseq object #####
    # Merge reads into Phyloseq object #
    #bac.unedited <- phyloseq::phyloseq(OTU.bac, TAX.bac, SAMP.bac, FASTA.bac, tree)
    ###IN ORDER TO KNIT THIS THE LINE ABOVE IS COMMENTED OUT; MULTIPLE FILES ARE COMPRESSED ###
    #save objects to reduce memory?
    #saveRDS(bac.unedited, file = "physeq_bac_full.rds")
    bac.unedited <- readRDS(file = "Bacteria_RDS/physeq_bac_full.rds")

# Decontaminate

    unique(bac.unedited@sam_data$Sample_Type)

    ## [1] "Spermosphere"       "Bulk_Soil"          "Positive_Control"  
    ## [4] "Extraction_Control" "PCR_Control"

    #bac.unedited@sam_data$Sample_or_Control <- ifelse(bac.unedited@sam_data$Sample_Type == "Positive_Control", "True Sample",
                                                      #ifelse(bac.unedited@sam_data$Sample_Type == "Negative_Control", "Control Sample",
                                                             #ifelse(bac.unedited@sam_data$Sample_Type == "PCR_Control", "Control Sample", "True Sample")))

    sample_data(bac.unedited)$is.neg <- sample_data(bac.unedited)$Sample_Type == "Extraction_Control" | sample_data(bac.unedited)$Sample_Type == "PCR_Control"
    contamdf.prev <- isContaminant(bac.unedited, method="prevalence", neg="is.neg", threshold = 0.1, normalize = TRUE)
    badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

    print(badTaxa) # 523 taxa

    ##   [1] "BOTU_114209" "BOTU_76094"  "BOTU_226201" "BOTU_180622" "BOTU_224189"
    ##   [6] "BOTU_214659" "BOTU_338722" "BOTU_106721" "BOTU_74675"  "BOTU_323940"
    ##  [11] "BOTU_287521" "BOTU_31225"  "BOTU_187079" "BOTU_229968" "BOTU_303340"
    ##  [16] "BOTU_92765"  "BOTU_256419" "BOTU_270593" "BOTU_100925" "BOTU_25821" 
    ##  [21] "BOTU_280843" "BOTU_147570" "BOTU_281732" "BOTU_102513" "BOTU_72786" 
    ##  [26] "BOTU_342012" "BOTU_9628"   "BOTU_92512"  "BOTU_221517" "BOTU_9762"  
    ##  [31] "BOTU_44628"  "BOTU_34583"  "BOTU_205217" "BOTU_271623" "BOTU_98696" 
    ##  [36] "BOTU_166810" "BOTU_82613"  "BOTU_340714" "BOTU_66240"  "BOTU_97150" 
    ##  [41] "BOTU_298633" "BOTU_115601" "BOTU_62722"  "BOTU_48485"  "BOTU_182960"
    ##  [46] "BOTU_32505"  "BOTU_32016"  "BOTU_180152" "BOTU_348268" "BOTU_133740"
    ##  [51] "BOTU_320198" "BOTU_103330" "BOTU_136450" "BOTU_321537" "BOTU_76674" 
    ##  [56] "BOTU_16873"  "BOTU_49352"  "BOTU_210420" "BOTU_85166"  "BOTU_227099"
    ##  [61] "BOTU_271450" "BOTU_246440" "BOTU_29901"  "BOTU_106771" "BOTU_57056" 
    ##  [66] "BOTU_31384"  "BOTU_13189"  "BOTU_306804" "BOTU_23018"  "BOTU_246518"
    ##  [71] "BOTU_128288" "BOTU_108341" "BOTU_81884"  "BOTU_154796" "BOTU_180793"
    ##  [76] "BOTU_339583" "BOTU_289289" "BOTU_259799" "BOTU_113152" "BOTU_121459"
    ##  [81] "BOTU_12215"  "BOTU_302221" "BOTU_28138"  "BOTU_128563" "BOTU_256169"
    ##  [86] "BOTU_233063" "BOTU_84579"  "BOTU_112090" "BOTU_217236" "BOTU_199430"
    ##  [91] "BOTU_229197" "BOTU_201835" "BOTU_7541"   "BOTU_46211"  "BOTU_22334" 
    ##  [96] "BOTU_168254" "BOTU_28406"  "BOTU_64347"  "BOTU_331609" "BOTU_91336" 
    ## [101] "BOTU_257369" "BOTU_48444"  "BOTU_89928"  "BOTU_272387" "BOTU_146992"
    ## [106] "BOTU_336283" "BOTU_89608"  "BOTU_178806" "BOTU_284242" "BOTU_155488"
    ## [111] "BOTU_290158" "BOTU_154360" "BOTU_82644"  "BOTU_36153"  "BOTU_349534"
    ## [116] "BOTU_91403"  "BOTU_258276" "BOTU_172108" "BOTU_59553"  "BOTU_194643"
    ## [121] "BOTU_256923" "BOTU_196944" "BOTU_229477" "BOTU_24095"  "BOTU_75671" 
    ## [126] "BOTU_182997" "BOTU_20062"  "BOTU_92389"  "BOTU_118781" "BOTU_35904" 
    ## [131] "BOTU_25739"  "BOTU_59047"  "BOTU_130514" "BOTU_320702" "BOTU_29685" 
    ## [136] "BOTU_167425" "BOTU_103936" "BOTU_105793" "BOTU_293375" "BOTU_13340" 
    ## [141] "BOTU_310302" "BOTU_67721"  "BOTU_60926"  "BOTU_200829" "BOTU_35111" 
    ## [146] "BOTU_23436"  "BOTU_55954"  "BOTU_164253" "BOTU_342190" "BOTU_315808"
    ## [151] "BOTU_88259"  "BOTU_74480"  "BOTU_323697" "BOTU_221313" "BOTU_181902"
    ## [156] "BOTU_192908" "BOTU_125462" "BOTU_266840" "BOTU_216974" "BOTU_46863" 
    ## [161] "BOTU_316257" "BOTU_293514" "BOTU_161297" "BOTU_187336" "BOTU_78611" 
    ## [166] "BOTU_155524" "BOTU_226935" "BOTU_258199" "BOTU_294033" "BOTU_322726"
    ## [171] "BOTU_315925" "BOTU_237200" "BOTU_304414" "BOTU_314103" "BOTU_11122" 
    ## [176] "BOTU_340864" "BOTU_299405" "BOTU_127039" "BOTU_317299" "BOTU_127736"
    ## [181] "BOTU_294511" "BOTU_254397" "BOTU_29815"  "BOTU_289605" "BOTU_33297" 
    ## [186] "BOTU_21665"  "BOTU_269930" "BOTU_62134"  "BOTU_216142" "BOTU_284832"
    ## [191] "BOTU_117710" "BOTU_63846"  "BOTU_48129"  "BOTU_222894" "BOTU_136452"
    ## [196] "BOTU_61215"  "BOTU_30523"  "BOTU_324977" "BOTU_229312" "BOTU_284179"
    ## [201] "BOTU_48058"  "BOTU_109754" "BOTU_40639"  "BOTU_133211" "BOTU_90043" 
    ## [206] "BOTU_196438" "BOTU_48681"  "BOTU_34353"  "BOTU_204504" "BOTU_43862" 
    ## [211] "BOTU_284621" "BOTU_133997" "BOTU_301853" "BOTU_339881" "BOTU_145899"
    ## [216] "BOTU_16310"  "BOTU_109874" "BOTU_18319"  "BOTU_348128" "BOTU_326664"
    ## [221] "BOTU_290957" "BOTU_282999" "BOTU_128718" "BOTU_25891"  "BOTU_42335" 
    ## [226] "BOTU_10152"  "BOTU_66342"  "BOTU_129618" "BOTU_115412" "BOTU_116752"
    ## [231] "BOTU_310614" "BOTU_282917" "BOTU_86066"  "BOTU_17229"  "BOTU_304265"
    ## [236] "BOTU_265577" "BOTU_157897" "BOTU_236694" "BOTU_90795"  "BOTU_168120"
    ## [241] "BOTU_157134" "BOTU_154339" "BOTU_213816" "BOTU_276793" "BOTU_68850" 
    ## [246] "BOTU_210219" "BOTU_272783" "BOTU_61222"  "BOTU_276726" "BOTU_193544"
    ## [251] "BOTU_22328"  "BOTU_156113" "BOTU_51972"  "BOTU_38652"  "BOTU_266896"
    ## [256] "BOTU_39577"  "BOTU_74836"  "BOTU_329566" "BOTU_47791"  "BOTU_73508" 
    ## [261] "BOTU_33434"  "BOTU_295988" "BOTU_261521" "BOTU_74875"  "BOTU_185465"
    ## [266] "BOTU_26789"  "BOTU_33443"  "BOTU_80369"  "BOTU_50177"  "BOTU_228701"
    ## [271] "BOTU_225516" "BOTU_267982" "BOTU_65153"  "BOTU_46809"  "BOTU_131436"
    ## [276] "BOTU_22069"  "BOTU_335881" "BOTU_33289"  "BOTU_101921" "BOTU_323120"
    ## [281] "BOTU_137929" "BOTU_83479"  "BOTU_176219" "BOTU_24568"  "BOTU_13807" 
    ## [286] "BOTU_161655" "BOTU_146834" "BOTU_315456" "BOTU_11000"  "BOTU_278200"
    ## [291] "BOTU_83287"  "BOTU_159214" "BOTU_123631" "BOTU_309347" "BOTU_298983"
    ## [296] "BOTU_162240" "BOTU_56316"  "BOTU_37771"  "BOTU_210530" "BOTU_33755" 
    ## [301] "BOTU_330934" "BOTU_293932" "BOTU_69783"  "BOTU_56163"  "BOTU_212035"
    ## [306] "BOTU_291902" "BOTU_54195"  "BOTU_40247"  "BOTU_312235" "BOTU_83777" 
    ## [311] "BOTU_57978"  "BOTU_201219" "BOTU_21486"  "BOTU_318518" "BOTU_10163" 
    ## [316] "BOTU_193428" "BOTU_54403"  "BOTU_62757"  "BOTU_281473" "BOTU_188644"
    ## [321] "BOTU_275859" "BOTU_65595"  "BOTU_284340" "BOTU_26870"  "BOTU_217904"
    ## [326] "BOTU_123581" "BOTU_141592" "BOTU_9692"   "BOTU_253423" "BOTU_200037"
    ## [331] "BOTU_13757"  "BOTU_44515"  "BOTU_139760" "BOTU_224362" "BOTU_189789"
    ## [336] "BOTU_254329" "BOTU_88284"  "BOTU_348306" "BOTU_145610" "BOTU_220986"
    ## [341] "BOTU_333018" "BOTU_79639"  "BOTU_37764"  "BOTU_219661" "BOTU_74221" 
    ## [346] "BOTU_239118" "BOTU_67814"  "BOTU_46712"  "BOTU_281806" "BOTU_96737" 
    ## [351] "BOTU_28788"  "BOTU_18530"  "BOTU_283881" "BOTU_345904" "BOTU_162605"
    ## [356] "BOTU_78123"  "BOTU_250742" "BOTU_311998" "BOTU_25384"  "BOTU_270791"
    ## [361] "BOTU_187539" "BOTU_91500"  "BOTU_25749"  "BOTU_203051" "BOTU_150636"
    ## [366] "BOTU_12332"  "BOTU_36565"  "BOTU_25826"  "BOTU_330736" "BOTU_254014"
    ## [371] "BOTU_35983"  "BOTU_140292" "BOTU_263473" "BOTU_346303" "BOTU_24203" 
    ## [376] "BOTU_270725" "BOTU_83500"  "BOTU_15866"  "BOTU_112296" "BOTU_345377"
    ## [381] "BOTU_173612" "BOTU_339864" "BOTU_55881"  "BOTU_19378"  "BOTU_34381" 
    ## [386] "BOTU_44921"  "BOTU_15079"  "BOTU_162272" "BOTU_187033" "BOTU_162539"
    ## [391] "BOTU_35345"  "BOTU_22621"  "BOTU_315472" "BOTU_173004" "BOTU_19191" 
    ## [396] "BOTU_86507"  "BOTU_13124"  "BOTU_56534"  "BOTU_344108" "BOTU_289621"
    ## [401] "BOTU_207721" "BOTU_73843"  "BOTU_231633" "BOTU_132107" "BOTU_260214"
    ## [406] "BOTU_234503" "BOTU_63062"  "BOTU_156754" "BOTU_272892" "BOTU_115128"
    ## [411] "BOTU_57938"  "BOTU_160852" "BOTU_39302"  "BOTU_150507" "BOTU_215315"
    ## [416] "BOTU_123673" "BOTU_54752"  "BOTU_211310" "BOTU_275567" "BOTU_310908"
    ## [421] "BOTU_58522"  "BOTU_49783"  "BOTU_273060" "BOTU_345406" "BOTU_288003"
    ## [426] "BOTU_287112" "BOTU_298824" "BOTU_147368" "BOTU_82436"  "BOTU_85178" 
    ## [431] "BOTU_204509" "BOTU_298146" "BOTU_169669" "BOTU_287565" "BOTU_162826"
    ## [436] "BOTU_190452" "BOTU_129500" "BOTU_206190" "BOTU_327575" "BOTU_254627"
    ## [441] "BOTU_107378" "BOTU_247837" "BOTU_115020" "BOTU_227502" "BOTU_186509"
    ## [446] "BOTU_121220" "BOTU_32489"  "BOTU_82578"  "BOTU_135910" "BOTU_226587"
    ## [451] "BOTU_286096" "BOTU_179556" "BOTU_16945"  "BOTU_321406" "BOTU_20687" 
    ## [456] "BOTU_42655"  "BOTU_265241" "BOTU_291401" "BOTU_106090" "BOTU_347158"
    ## [461] "BOTU_126438" "BOTU_205916" "BOTU_77429"  "BOTU_22569"  "BOTU_74568" 
    ## [466] "BOTU_140984" "BOTU_242317" "BOTU_302519" "BOTU_250539" "BOTU_140097"
    ## [471] "BOTU_49392"  "BOTU_125567" "BOTU_343288" "BOTU_76454"  "BOTU_341464"
    ## [476] "BOTU_83470"  "BOTU_28725"  "BOTU_50442"  "BOTU_35460"  "BOTU_322340"
    ## [481] "BOTU_280912" "BOTU_80883"  "BOTU_104187" "BOTU_231953" "BOTU_164961"
    ## [486] "BOTU_257889" "BOTU_319086" "BOTU_113161" "BOTU_56076"  "BOTU_50701" 
    ## [491] "BOTU_110686" "BOTU_136369" "BOTU_69991"  "BOTU_130768" "BOTU_239491"
    ## [496] "BOTU_250673" "BOTU_226849" "BOTU_28646"  "BOTU_319467" "BOTU_193660"
    ## [501] "BOTU_129520" "BOTU_277029" "BOTU_311290" "BOTU_278389" "BOTU_143105"
    ## [506] "BOTU_204641" "BOTU_324659" "BOTU_173939" "BOTU_312189" "BOTU_45018" 
    ## [511] "BOTU_39831"  "BOTU_265939" "BOTU_190468" "BOTU_132117" "BOTU_131314"
    ## [516] "BOTU_325973" "BOTU_42935"  "BOTU_148436" "BOTU_82986"  "BOTU_173846"
    ## [521] "BOTU_23402"  "BOTU_320289" "BOTU_121069"

    ps.pa <- transform_sample_counts(bac.unedited, function(abund) 1*(abund>0))
    ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg == "TRUE", ps.pa)
    ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == "FALSE", ps.pa)

    view(ps.pa.pos@sam_data)

    # Make data.frame of prevalence in positive and negative samples
    df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                        contaminant=contamdf.prev$contaminant)
    decontaminate.bac <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
      geom_point() +
      xlab("Prevalence (Negative Controls)") + 
      ylab("Prevalence (True Samples)") + 
      scale_color_manual(values = cbbPalette)+ 
      ggtitle("Prokaryote") +
      theme_classic()
    goodTaxa <- setdiff(taxa_names(bac.unedited), badTaxa)
    str(goodTaxa) #124,214 good taxa

    ##  chr [1:349319] "BOTU_271591" "BOTU_167044" "BOTU_238267" "BOTU_323756" ...

    bac_filt <- prune_taxa(goodTaxa, bac.unedited)

# Taxonomy filtering

    # remove OTUs that are mitochondria, chloroplast, or unidentified at the kingdom level 
    bac_no_chloro <- bac_filt %>% 
      phyloseq::subset_taxa(Order != "Chloroplast") %>%
      phyloseq::subset_taxa(Family != "Mitochondria") %>%
      phyloseq::subset_taxa(Kingdom != "unidentified") %>%
      core(detection = 2, prevalence = 0.20)

    # Number of reads to Chloroplast, Mitochondria, and unidentified
    chloro.mito.reads <- sum(sample_sums(bac.unedited))-sum(sample_sums(bac_no_chloro)) #17660 Reads

    # Percent reads that were chloroplast, mito, or unidentified
    (chloro.mito.reads/sum(sample_sums(bac.unedited)))*100

    ## [1] 37.374

    #20.1 percent

# Mock Community analysis

    # positive controls
    bac_mock <- bac_no_chloro %>% 
      subset_samples(Sample_Type == "Positive_Control") #%>%
      #phyloseq::filter_taxa(function(x) sum(x) > 1500, TRUE) # filter OTUs to have more than 10 read in mock samples

    mock2 <- microbiome::transform(bac_mock, "compositional") # relative abundance transform

    sequenced.mock.bac <- mock2 %>%
      psmelt() %>% 
      ggplot(aes(Sample, Abundance, fill = Genus)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      #scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "", y = "Relative abundance (%)",
           title = "Prokaryote") + 
      theme(axis.text.x = element_text(angle=45, hjust=1)) 

    ## Warning in psmelt(.): The rank names: 
    ## OTU
    ##  have been renamed to: 
    ## taxa_OTU
    ## to avoid conflicts with special phyloseq plot attribute names.

    sequenced.mock.bac

![](LL_Prokaryote_Publication_Preprocessing_files/figure-markdown_strict/Mock%20Community%20Analysis-1.png)

    # Adding in theoretical distribution - the last two are fungi and are not expected to be amplified with 16S
    Label <- c("Pseudomonas aeruginosa", 
               "Escherichia coli",
               "Salmonella enterica", 
               "Lactobacillus fermentum", 
               "Enterococcus faecalis", 
               "Staphylococcus aureus", 
               "Listeria monocytogenes", 
               "Bacillus subtilis")

    # theoretical species composition in the mock community
    Abundance <- c(rep(0.125, 8))

    th.mock <- data.frame(Label, Abundance)
    th.mock$Sample <- "Theoretical"

    th.mock$Label <- factor(th.mock$Label, levels = c("Lactobacillus fermentum", 
                                                      "Staphylococcus aureus", 
                                                      "Bacillus subtilis",
                                                      "Escherichia coli",
                                                      "Listeria monocytogenes",
                                                      "Enterococcus faecalis",
                                                      "Salmonella enterica",
                                                      "Pseudomonas aeruginosa"))


    theory.mock <- ggplot(th.mock, aes(Sample, Abundance, fill = Label)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      scale_fill_manual(values= c(cbbPalette[[1]], 
                                  cbbPalette[[2]], 
                                  cbbPalette[[3]], 
                                  cbbPalette[[4]], 
                                  cbbPalette[[5]],
                                  cbbPalette[[6]],
                                  cbbPalette[[8]],
                                  "violet", "pink", "grey", "black", "blue")) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "", y = "Relative abundance (%)",
           title = "Theoretical composition") + 
      theme(axis.text.x = element_text(angle=45, hjust=1),
            legend.text = element_text(face = "italic"),
            legend.title = element_blank())

    # I think maybe the theoretical mock community can also be mentioned in the figure legend. 

    mock.composition <- mock2 %>%
      psmelt() %>%
      dplyr::group_by(OTU) %>%
      dplyr::summarise(MeanRelAbund = mean(Abundance)) %>%
      arrange(-MeanRelAbund)

    ## Warning in psmelt(.): The rank names: 
    ## OTU
    ##  have been renamed to: 
    ## taxa_OTU
    ## to avoid conflicts with special phyloseq plot attribute names.

    # these 8 OTUs made up 99.3% of the mock composition. These OTUs also match the 8 supposed to be in the mock
    sum(mock.composition[1:8,]$MeanRelAbund)

    ## [1] 0.9496684

# Data filtering

    # remove samples with less than 5000 reads
    bac_sperm <- bac_no_chloro %>% 
      subset_samples(Sample_Type %in% c("Bulk_Soil", "Spermosphere")) #%>%
      #prune_samples(sample_sums(.) > 5000, .) %>% # remove samples below 5000 reads
      #phyloseq::filter_taxa(function(x) sum(x) > 1, TRUE) # remove taxa with less than 2 reads
    sample.sums <- data.frame(sample_sums(bac_no_chloro))
    #One sample was removed due to low sequencing coverage

    ###### RDS of Non-normalized Prokaryote data ######
    # Save an object to a file
    saveRDS(bac_sperm, file = "Bacteria_RDS/Spermosphere_bac_clean.rds")

# Read Depth Analysis

    # Restore the object
    bac_sperm <- readRDS(file = "Bacteria_RDS/Spermosphere_bac_clean.rds")

    ###### READS PER SAMPLE ######
    sample.sums <- data.frame(sample_sums(bac_sperm))

    view(sample.sums)

    read.dist.bac <- ggplot(sample.sums, aes(x = sample_sums.bac_sperm.)) +
      geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
      theme_classic() +
      xlab("Read Depth") + 
      ggtitle("Prokaryote")

    sum(sample_sums(bac_sperm)) # total reads = 38960707

    ## [1] 30744771

    median(sample_sums(bac_sperm)) # 100489.5

    ## [1] 76247.5

    str(bac_sperm@sam_data) #384 Samples

    ## 'data.frame':    384 obs. of  33 variables:
    ## Formal class 'sample_data' [package "phyloseq"] with 4 slots
    ##   ..@ .Data    :List of 33
    ##   .. ..$ : int  201 202 203 204 205 206 207 208 209 210 ...
    ##   .. ..$ : chr  "A11" "A12" "A13" "A14" ...
    ##   .. ..$ : chr  "0.2157" "0.1736" "0.1892" "0.1702" ...
    ##   .. ..$ : chr  "0.2157" "0.1736" "0.1892" "0.1702" ...
    ##   .. ..$ : chr  "0" "0" "0" "0" ...
    ##   .. ..$ : chr  "0" "0" "0" "0" ...
    ##   .. ..$ : chr  "" "" "" "" ...
    ##   .. ..$ : chr  "0" "0" "0" "0" ...
    ##   .. ..$ : chr  "4/3/2024" "4/3/2024" "4/3/2024" "4/3/2024" ...
    ##   .. ..$ : chr  "4/3/2024" "4/3/2024" "4/3/2024" "4/3/2024" ...
    ##   .. ..$ : chr  "1/31/2025" "1/31/2025" "1/31/2025" "1/31/2025" ...
    ##   .. ..$ : chr  "3/12/2025" "3/12/2025" "3/12/2025" "3/12/2025" ...
    ##   .. ..$ : chr  "bc43" "bc44" "bc45" "bc46" ...
    ##   .. ..$ : chr  "0.121" "0.14" "0.193" "0.219" ...
    ##   .. ..$ : chr  "too low" "too low" "too low" "too low" ...
    ##   .. ..$ : chr  "Autoclaved_Potting_Soil" "Autoclaved_Potting_Soil" "Autoclaved_Potting_Soil" "Autoclaved_Potting_Soil" ...
    ##   .. ..$ : chr  "Control" "Control" "Control" "Control" ...
    ##   .. ..$ : chr  "NA" "NA" "NA" "NA" ...
    ##   .. ..$ : chr  "NA" "NA" "NA" "NA" ...
    ##   .. ..$ : chr  "NA" "NA" "NA" "NA" ...
    ##   .. ..$ : chr  "Spermosphere" "Spermosphere" "Spermosphere" "Spermosphere" ...
    ##   .. ..$ : chr  "Control" "Control" "Control" "Control" ...
    ##   .. ..$ : chr  "Control" "Control" "Control" "Control" ...
    ##   .. ..$ : chr  "Control" "Control" "Control" "Control" ...
    ##   .. ..$ : chr  "" "" "" "" ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  NA NA NA NA NA NA ...
    ##   .. ..$ : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##   ..@ names    : chr  "Tube_Number" "Sample_ID" "Presowing_Weight" "Postsowing_Weight" ...
    ##   ..@ row.names: chr  "A11" "A12" "A13" "A14" ...
    ##   ..@ .S3Class : chr "data.frame"

    str(bac_sperm@otu_table) #123518 OTU

    ## Formal class 'otu_table' [package "phyloseq"] with 2 slots
    ##   ..@ .Data        : int [1:4947, 1:384] 0 0 1 0 0 0 0 0 0 0 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. ..$ : chr [1:4947] "BOTU_1424" "BOTU_1762" "BOTU_620" "BOTU_3916" ...
    ##   .. .. ..$ : chr [1:384] "A11" "A12" "A13" "A14" ...
    ##   ..@ taxa_are_rows: logi TRUE
    ##   ..$ dim     : int [1:2] 4947 384
    ##   ..$ dimnames:List of 2
    ##   .. ..$ : chr [1:4947] "BOTU_1424" "BOTU_1762" "BOTU_620" "BOTU_3916" ...
    ##   .. ..$ : chr [1:384] "A11" "A12" "A13" "A14" ...

# Rarefaction Analysis

    ###### Rarefaction analysis #####
    sam.data <- data.frame(bac_sperm@sam_data)
    bOTU.table <- otu_table(bac_sperm) %>%
      as.data.frame() %>%
      as.matrix()

    raremax <- min(rowSums(t(bOTU.table)))
    raremax2 <- min(rowSums(bOTU.table))

    rare.fun <- rarecurve(t(bOTU.table), step = 1000, sample = raremax, tidy = T)

    bac.rare.curve.extract2 <- left_join(sam.data, rare.fun, by = c("Sample_ID" = "Site"))

    bac.rare <- ggplot(bac.rare.curve.extract2, aes(x = Sample, y = Species, group = Sample_ID, color = Sample_Type)) + 
      #geom_point() +
      geom_line() + 
      xlab("Reads") + 
      ylab("Number of OTUs") +
      ggtitle("Prokaryote") +
      theme_classic() + 
      geom_vline(xintercept = median(sample_sums(bac_sperm)), linetype = "dashed") +
      scale_color_manual(values = cbbPalette)

    bac.rare2 <- ggplot(bac.rare.curve.extract2, aes(x = Sample, y = Species, group = Sample_ID, color = Location_EC)) + 
      #geom_point() +
      geom_line() + 
      xlab("Reads") + 
      ylab("Number of OTUs") +
      ggtitle("Prokaryote") +
      theme_classic() + 
      geom_vline(xintercept = median(sample_sums(bac_sperm)), linetype = "dashed") +
      scale_color_manual(values = cbbPalette)

    ggarrange(bac.rare, 
              read.dist.bac, 
              sequenced.mock.bac, 
              decontaminate.bac, nrow = 2, ncol = 2, labels = "auto")

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](LL_Prokaryote_Publication_Preprocessing_files/figure-markdown_strict/Rarefaction%20Analysis-1.png)

    ggarrange(bac.rare, bac.rare2, nrow = 1, ncol =2)

![](LL_Prokaryote_Publication_Preprocessing_files/figure-markdown_strict/Rarefaction%20Analysis-2.png)

# Metagenome CSS normalization

    MGS <- phyloseq_to_metagenomeSeq(bac_sperm) #converts to metagenomeseq format
    p <- metagenomeSeq::cumNormStatFast(MGS)

    ## Default value being used.

    MGS <- metagenomeSeq::cumNorm(MGS, p =p)
    metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample

    ##   A11   A12   A13   A14   A15   A16   A21   A22   A23   A24   A25   A26   A31 
    ##   874   385   489   260   680   452   916   374   904  2160  2316  1024   906 
    ##   A32   A33   A34   A35   A36   A41   A42   A43   A44   A45   A46  BA11  BA12 
    ##  2073  2296  2005  1043   945   322   960  2200   952   944   912   471   435 
    ##  BA13  BA14  BA15  BA16  BA21  BA22  BA23  BA24  BA25  BA26  BA31  BA32  BA33 
    ##   398   618   441   918   375   409   533   615   243   328   286   802   428 
    ##  BA34  BA35  BA36  BA41  BA42  BA43  BA44  BA45  BA46  BE11  BE12  BE13  BE14 
    ##   619   445   380   261   313   537   491   276   392  8604  8587  7664  7526 
    ##  BE15  BE16  BE21  BE22  BE23  BE24  BE25  BE26  BE31  BE32  BE33  BE34  BE35 
    ## 11313  9831  8528  9582  7508  8747  8628  8509  8718  6625  7327  9705  9519 
    ##  BE36  BE41  BE42  BE43  BE44  BE45  BE46  BM11  BM12  BM13  BM14  BM15  BM16 
    ##  9628  7449  9594  9587  5059  7566  6311  6300  5965  5229  4946  5799  5756 
    ##  BM21  BM22  BM23  BM24  BM25  BM26  BM31  BM32  BM33  BM34  BM35  BM36  BM41 
    ##  4916  6298  5001  5804  3948  6415  6191  4803  6057  5037  5801  5012  5836 
    ##  BM42  BM43  BM44  BM45  BM46  BN11  BN12  BN13  BN14  BN15  BN16  BN21  BN22 
    ##  4934  4998  5786  3909  4906  1548  1754  1460  1505  1374  2067  1710  1877 
    ##  BN23  BN24  BN25  BN26  BN31  BN32  BN33  BN34  BN35  BN36  BN41  BN42  BN43 
    ##  1837  1544  1658  1657  1663  2346  1808  1824  1249  1237  1725  1882  2209 
    ##  BN44  BN45  BN46  BT11  BT12  BT13  BT14  BT15  BT16  BT21  BT22  BT23  BT24 
    ##  1289  1933  1448  5215  4362  5222  4755  4612  4351  7101  9256  5471  5586 
    ##  BT25  BT26  BT31  BT32  BT33  BT34  BT35  BT36  BT41  BT42  BT43  BT44  BT45 
    ##  6462  5743  4436  7238  6065  6276  7306  6980  4920  3923  3874  4796  3699 
    ##  BT46  BW11  BW12  BW13  BW14  BW15  BW16  BW21  BW22  BW23  BW24  BW25  BW26 
    ##  3789  7250  7072  7048  6813  5952  6036  6692  7020  6055  5855  9195  6835 
    ##  BW31  BW32  BW33  BW34  BW35  BW36  BW41  BW42  BW43  BW44  BW45  BW46   E11 
    ##  5911  8083  6782  7028  6627  8927  7350  7104  6843  9023  6078  8170  7853 
    ##   E12   E13   E14   E15   E16   E21   E22   E23   E24   E25   E26   E31   E32 
    ##  7537  6546  6698  6851  8447  6306  5363  5346  6278  6417  4317  6573  8048 
    ##   E33   E34   E35   E36   E41   E42   E43   E44   E45   E46   M11   M12   M13 
    ##  6704  6646  6642  7798  6662  5604  6585  9273  7927  4420  6758  4917  6338 
    ##   M14   M15   M16   M21   M22   M23   M24   M25   M26   M31   M32   M33   M34 
    ##  3982  5131  4999  3262  4108  6521  4083  5241  4293  2816  3822  5099  6302 
    ##   M35   M36   M41   M42   M43   M44   M45   M46   N11   N12   N13   N14   N15 
    ##  4199  1938  5086  2972  4097  2979  4986  4081  2278  1370  1347  2437  1733 
    ##   N16   N21   N22   N23   N24   N25   N26   N31   N32   N33   N34   N35   N36 
    ##  1371  1717  1367  1514  1913  1717  1849  1729  1322  1528  1401  1833  1976 
    ##   N41   N42   N43   N44   N45   N46   P11   P12   P13   P14   P15   P16   P21 
    ##  1814  1237  1784  1309  1234  1320  5624  4550  5583  4744  3777  6394  5308 
    ##   P22   P23   P24   P25   P26   P31   P32   P33   P34   P35   P36   P41   P42 
    ##  4761  5747  4566  4083  4850  3813  5843  6764  6120  3922  3686  3820  2703 
    ##   P43   P44   P45   P46  PB11  PB12  PB13  PB14  PB15  PB16  PB21  PB22  PB23 
    ##  5023  3484  2755  3099  5467  5641  4449  4529  5408  5453  4468  3679  4316 
    ##  PB24  PB25  PB26  PB31  PB32  PB33  PB34  PB35  PB36  PB41  PB42  PB43  PB44 
    ##  5728  4496  4680  4411  4374  4382  4666  5612  3614  4471  5488  4823  4800 
    ##  PB45  PB46   R11   R12   R13   R14   R15   R16   R21   R22   R23   R24   R25 
    ##  3752  4851  3486  3420  3736  4530  3539  3487  2212  3582  3257  3365  3261 
    ##   R26   R31   R32   R33   R34   R35   R36   R41   R42   R43   R44   R45   R46 
    ##  2712  2565  3302  2553  3724  2553  2526  3791  2560  2430  3488  2700  2438 
    ##  RB11  RB12  RB13  RB14  RB15  RB16  RB21  RB22  RB23  RB24  RB25  RB26  RB31 
    ##  3734  4783  3502  4477  3485  3675  4595  3833  2903  1587  2760  3491  4369 
    ##  RB32  RB33  RB34  RB35  RB36  RB41  RB42  RB43  RB44  RB45  RB46   T11   T12 
    ##  5660  3658  4686  4227  4290  5396  5038  5500  4227  4378  5612  5395  6553 
    ##   T13   T14   T15   T16   T21   T22   T23   T24   T25   T26   T31   T32   T33 
    ##  3707  6987  2782  3778  3759  3544  3640  4812  2651  4706  6216  4683  3802 
    ##   T34   T35   T36   T41   T42   T43   T44   T45   T46   W11   W12   W13   W14 
    ##  4651  4600  3703  4480  4919  3673  3593  4138  4547  8144  7199  8590  7058 
    ##   W15   W16   W21   W22   W23   W24   W25   W26   W31   W32   W33   W34   W35 
    ##  6477  6144  8283  8183  7459  8247  7187  8413  6087  4056  8164  6254  6262 
    ##   W36   W41   W42   W43   W44   W45   W46 
    ##  6151  7158  6382  6196  7974  7237  4221

    norm.bac <- metagenomeSeq::MRcounts(MGS, norm = T) 
    norm.bac.OTU <- phyloseq::otu_table(norm.bac, taxa_are_rows = TRUE) #exports the new otu table

    #bac.css.norm <- phyloseq::phyloseq(norm.bac.OTU, FASTA.bac, SAMP.bac, TAX.bac, tree) #new otu table phyloseq object
    # ABOVE LINE IS COMMENTED OUT BECAUSE THE FILE IS ZIPPED IN A FILE #
    #saveRDS(bac.css.norm, file = "Bacteria_RDS/Spermosphere_CSS.rds")
    # ABOVE LINE IS COMMENTED OUT BECAUSE THE FILE IS ZIPPED IN A FILE #

    # Restore the object
    bac.css.norm <- readRDS(file = "Bacteria_RDS/Spermosphere_CSS.rds")



    bs.rarefied <- rarefy_even_depth(bac_sperm, rngseed=12345, sample.size=0.9*min(sample_sums(bac_sperm)), replace=F)

    ## `set.seed(12345)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(12345); .Random.seed` for the full vector

    ## ...

    bs.rarefied

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 4947 taxa and 384 samples ]
    ## sample_data() Sample Data:       [ 384 samples by 33 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 4947 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 4947 tips and 4945 internal nodes ]
    ## refseq()      DNAStringSet:      [ 4947 reference sequences ]

    saveRDS(bs.rarefied, file = "Spermosphere_Rare.rds")
