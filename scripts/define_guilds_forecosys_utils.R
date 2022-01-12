read_coverage <- function(genome2representative.file, coveragelong.file, metadata.file) {
  # ecosysguilds = read.table(file.path("/Users/ukaraoz/Work/Isogenie/IsogenieGenomes/permafrost-modelling/results", "IsogenieGenomes.ecosysguilds.xls"), sep = "\t")
  # # 1,529 Stordalen MAGs dereplicated into 647 genomes
  # genomes_ecosysguilds = setdiff(ecosysguilds[,1], c("genome","")) %>% tibble::as_tibble() %>% dplyr::pull(value)
  # genome2representativeGenomes = read.table(genome2representative.file, header = T, check.names = F, sep = "\t") %>% 
  #   tibble::as_tibble() %>%
  #   dplyr::pull(c("Genome")) %>% unique()
  # genome2representativeClusterReps = read.table(genome2representative.file, header = T, check.names = F, sep = "\t") %>% 
  #     tibble::as_tibble() %>%
  #     dplyr::pull(c("97% ANI cluster representative")) %>% unique()

  # 647 representatives
  genome2representative = read.table(genome2representative.file, header = T, check.names = F, sep = "\t") %>% 
    tibble::as_tibble() %>%
    dplyr::select(c("Genome", "97% ANI cluster representative", "Clade name", "Named subclade"))
  
  # read coverage
  coveragelong = read.table(coveragelong.file, header = T, sep = "\t") %>% tibble::as_tibble()
  sample2metadata = read.table(metadata.file, header = T, sep = "\t") %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(samplelong = paste0("month=", month, "_",
                                      "year=", year, "_", 
                                      "habitat=", habitat, "_",
                                      "core=", core, "_",
                                      "depth=", depth, "_",
                                      "id=", sample)) %>%
    dplyr::select(c("sample", "samplelong"))
  coveragelong = coveragelong %>% dplyr::left_join(sample2metadata, by = c("sample" = "sample"))
  coveragewide = coveragelong %>% 
    dplyr::select(c("genome", "samplelong", "coverage")) %>% 
    tidyr::pivot_wider(names_from = samplelong, values_from = coverage) %>% 
    # join coverage with cluster representative info
    dplyr::left_join(genome2representative, by = c("genome" = "Genome")) %>% 
    dplyr::select(-c("genome", "Clade name", "Named subclade")) %>%
    #dplyr::select(c("genome", "97% ANI cluster representative", "Clade name", "Named subclade"), everything())
    dplyr::select(c("97% ANI cluster representative"), everything()) %>%
    dplyr::mutate(`ANI97%clusterrepresentative` = `97% ANI cluster representative`)

  relabundwide = coveragelong %>% 
    dplyr::select(c("genome", "samplelong", "relative_abundance")) %>% 
    tidyr::pivot_wider(names_from = samplelong, values_from = relative_abundance) %>% 
    # join coverage with cluster representative info
    dplyr::left_join(genome2representative, by = c("genome" = "Genome")) %>% 
    dplyr::select(-c("genome", "Clade name", "Named subclade")) %>%
    #dplyr::select(c("genome", "97% ANI cluster representative", "Clade name", "Named subclade"), everything())
    dplyr::select(c("97% ANI cluster representative"), everything()) %>%
    dplyr::mutate(`ANI97%clusterrepresentative` = `97% ANI cluster representative`)

  relabundofrecovered = coveragelong %>% 
    dplyr::select(c("genome", "samplelong", "relative_abundance_of_recovered")) %>% 
    tidyr::pivot_wider(names_from = samplelong, values_from = relative_abundance_of_recovered) %>% 
    # join coverage with cluster representative info
    dplyr::left_join(genome2representative, by = c("genome" = "Genome")) %>% 
    dplyr::select(-c("genome", "Clade name", "Named subclade")) %>%
    #dplyr::select(c("genome", "97% ANI cluster representative", "Clade name", "Named subclade"), everything())
    dplyr::select(c("97% ANI cluster representative"), everything()) %>%
    dplyr::mutate(`ANI97%clusterrepresentative` = `97% ANI cluster representative`)

  result = list()
  result$coverage = coveragewide
  result$relabund = relabundwide
  result$relabundofrecovered = relabundofrecovered
  result
}

build_metadata <- function(taxonomy.file, metadata.file) {
  taxonomy = read.table(taxonomy.file, header = T, sep = "\t") %>% tibble::as_tibble() %>%
    dplyr::mutate(sample = stringr::str_replace(`genome`, "(.*)\\.(.*)\\.(.*)", "\\2"))
  sample2metadata = read.table(metadata.file, header = T, sep = "\t") %>% tibble::as_tibble()
  sample2taxonomymetadata = taxonomy %>%
    dplyr::inner_join(sample2metadata, by = c("sample" = "sample"))
  sample2taxonomymetadata
}

select_traitsrulesgenes <- function() {
  denitrification_traits = 
    c("Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:N (denitrification):denitrification to nitrite",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:N (denitrification):denitrification to nitric oxide",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:N (denitrification):denitrification to nitrous oxide",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:N (denitrification):denitrification to nitrogen")
  NH4_NO2_oxidation_traits = 
    c("Resource Use:Chemotrophy:chemolithoautotrophy:aerobic ammonia oxidation",
      "Resource Use:Chemotrophy:chemolithoautotrophy:nitrite oxidation")
  nitrogenfixation_traits = c("Resource Acquisition:Substrate assimilation:N compounds:nitrogen fixation")
  methanemethanoloxidation_traits = 
    c("Resource Acquisition:Substrate assimilation:C1 compounds:methanotrophy:methane oxidation pathway",
      "Resource Acquisition:Substrate assimilation:C1 compounds:methanotrophy:methanol oxidation pathway")
  acetoclasticmethanogenesis_traits = c("Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:acetoclastic methanogenesis")
  hydrogenotrophicmethanogenesis_traits = c("Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:hydrogenotrophic methanogenesis")
  methylotrophicmethanogenesis_traits = 
    c("Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:methylotrophic methanogenesis:from methanethiol",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:methylotrophic methanogenesis:from methanol",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:methylotrophic methanogenesis:from TMAO",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:methylotrophic methanogenesis:from monomethylamine",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:methylotrophic methanogenesis:from dimethylamine",
      "Resource Use:Chemotrophy:chemolithoheterotrophy:anaerobic respiration:methylotrophic methanogenesis:from trimethylamine")
  fermentation_traits = 
    c("Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of pyruvate to acetone",
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of pyruvate to lactate and ethanol", 
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of pyruvate to alcohols", 
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of pyruvate to SCFAs",
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of pyruvate to butylene glycol", 
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of pyruvate to isobutanol", 
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of acetylene to acetate", 
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:fermentation of acetylene to acetate and ethanol")
  heterotrophic_traits = 
    c("Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:ED pathway", 
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:EMP pathway",
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:fructose degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:fucose degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:maltose degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:galactose degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:mannose degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:carbohydrate degradation:trehalose degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:sugar acid degradation:galacturonate degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:sugar acid degradation:glucuronate degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:alcohol degradation:glycerol degradation",
    "Resource Acquisition:Substrate degradation:simple compound degradation:phenol degradation:resorcinol degradation")
  ETC_traits = c("Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex I",
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex II", 
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III",
    "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex IV")
  select_traits = 
    c(denitrification_traits, 
      NH4_NO2_oxidation_traits, 
      nitrogenfixation_traits, 
      methanemethanoloxidation_traits,
      hydrogenotrophicmethanogenesis_traits,
      acetoclasticmethanogenesis_traits,
      methylotrophicmethanogenesis_traits,
      fermentation_traits,
      heterotrophic_traits, 
      ETC_traits)

  # rules
  denitrification_rules = c("nitrite->nitricoxide_partial", "nitricoxide->nitrousoxide_partial", "nitrousoxide->dinitrogen")  
  NH4_NO2_oxidation_rules = c("pmoamoABC_or", "hydroxylamine->nitrite", "nitrite->nitrate_partial")
  nitrogenfixation_rules = c("dinitrogen->ammonia_partial")
  methanemethanoloxidation_rules = c("CH4->CH3OH", "CH3OH->CH2O")
  acetoclasticmethanogenesis_rules = c("acetyl-P<->acetate", "acetyl-P->acetyl-CoA", "acetyl-CoA+methylH4MPT->methylH4SPT+CO_partial", "methylCoM+CoB->CoM+CH4_or")
  hydrogenotrophicmethanogenesis_rules = c("MF+CO2->formylMF_or", "formylMF+H4MPT->MF+formylH4MPT", "formylH4MPT->methenylH4MPT", "methenylH4MPT->methyleneH4MPT", "methyleneH4MPT->methyl-H4SPT", "methylCoM+CoB->CoM+CH4_or")
  methylotrophicmethanogenesis_rules = c("methanethiol->methylCoM", "methanol->methylCoM_or", "TMAO->trimethylamine", "monomethylamine->methylCoM-direct_or", "dimethylamine->methylCoM", "trimethylamine->methylCoM", "methylCoM+CoB->CoM+CH4_or")
  fermentation_rules = c("acetoacetyl-CoA->acetone", "pyruvate->lactate and ethanol", "pyruvate->alcohols", "pyruvate->SCFAs", "pyruvate->butyleneglycol", "pyruvate->butanol", "acetylene->acetate", "acetylene->acetate+ethanol")
  heterotrophic_rules = c("glycolysis-ED_full", "glycolysis-EMP_full", "fructose degradation", "fucose degradation", "maltose degradation", "galactose degradation", "mannose degradation", "trehalose degradation", "galacturonate degradation", "glucuronate degradation", "glycerol degradation", "resorcinol degradation")
  ETC_rules = c("nuo_and", "succinate->fumarate_key", "cytochromecreductase_or", "cytochromeoubiquinol_and", "cytochromec_and", "cytochromeaa3600menaquinol_and", "cytochromeccbb3oxidase_and")

  select_rules = c(denitrification_rules, NH4_NO2_oxidation_rules, nitrogenfixation_rules,
                   nitrogenfixation_rules, methanemethanoloxidation_rules, 
                   acetoclasticmethanogenesis_rules, hydrogenotrophicmethanogenesis_rules,
                   methylotrophicmethanogenesis_rules, fermentation_rules, heterotrophic_rules,
                   ETC_rules)

  denitrification_genes = c("narG", "narH", "narI_dsrM", "napA", "napB", "nirS", "nirK", "norB", "norC", "norV", "norW", "nosZ")
  NH4_NO2_oxidation_genes = c("pmoA-amoA", "pmoB-amoB", "pmoC-amoC", "hao", "narG", "narH")
  nitrogenfixation_genes = c("nifD", "nifK", "nifH", "anfG", "vnfD", "vnfG", "vnfK", "vnfH")
  methaneoxidation_genes = c("mmoX", "mmoY", "mmoZ", "mmoC", "mmoD", "pmoA-amoA", "pmoB-amoB", "pmoC-amoC")
  methanoloxidation_genes = c("mxaF", "mxaJ", "mxaG", "mxaI", "mxaA", "mxaC", "mxaD", "mxaK", "mxaL", "xoxF")
  methanogenesislowerbranch_genes = c("tmtrA", "tmtrB", "tmtrC", "tmtrD", "tmtrE", "tmtrF", "tmtrG", "tmtrH", "mcrA", "mcrB", "mcrG", "mcr2", "mcrC", "mcrD")
  methanogenesisacetoclastic_genes = c("ackA", "pta1", "pta2", "cdhC", "acsC", "acsD")
  methanogenesishydrogenotrophic_genes = c("fwdA", "fwdB", "fwdC", "fwdD", "fwdF", "fwdH", "fwdG", "fwdE", "ftr", "mtdch", "hmd", "mer")
  methanogenesismethylotrophic_genes = c("mmtsA", "mtaB", "mtaC", "torA", "torZ", "torC", "torY", "mtbA", "mtmB", "mtmC", "mtbA", "mtbB", "mtbC", "mtbA", "mttB", "mttC")
  ETCI_genes = c("nuoA", "nuoB", "nuoBCD", "nuoC", "nuoCD", "nuoD", "nuoE", "nuoF", "nuoG", "nuoH", "nuoI", "nuoJ", "nuoK", "nuoL", "nuoLM", "nuoM", "nuoN")
  ETCII_genes = c("sdhA", "sdhB", "sdhC", "sdhD", "frdA", "frdB", "frdC", "frdD", "sdhD_arch1", "sdhD_arch2")
  ETCIII_genes = c("qcrA", "fbcH", "petB_2", "petB_1", "qcrB")
  ETCIV_genes = c("cyoA", "cyoB", "cyoC", "coxA", "coxB", "coxC", "qoxA", "qoxB", "qoxC", "ccoN", "ccoO", "ccoP")
  
  select_genes = c(denitrification_genes, NH4_NO2_oxidation_genes, nitrogenfixation_genes, 
                   methaneoxidation_genes, methanoloxidation_genes,
                   methanogenesislowerbranch_genes, methanogenesisacetoclastic_genes, 
                   methanogenesishydrogenotrophic_genes, methanogenesismethylotrophic_genes,
                   ETCI_genes, ETCII_genes, ETCIII_genes, ETCIV_genes)
  results = list(denitrification_traits = denitrification_traits,
                 NH4_NO2_oxidation_traits = NH4_NO2_oxidation_traits,
                 nitrogenfixation_traits = nitrogenfixation_traits,
                 methanemethanoloxidation_traits = methanemethanoloxidation_traits,
                 hydrogenotrophicmethanogenesis_traits = hydrogenotrophicmethanogenesis_traits,
                 acetoclasticmethanogenesis_traits = acetoclasticmethanogenesis_traits,
                 methylotrophicmethanogenesis_traits = methylotrophicmethanogenesis_traits,
                 fermentation_traits = fermentation_traits,
                 heterotrophic_traits = heterotrophic_traits,
                 ETC_traits = ETC_traits,
                 select_traits = select_traits, 
                 select_rules = select_rules, 
                 select_genes = select_genes)
  results
}