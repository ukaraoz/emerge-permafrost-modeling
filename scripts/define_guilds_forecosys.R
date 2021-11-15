library(dplyr)
base = "/Users/ukaraoz/Work/Isogenie/IsogenieGenomes"
source("/Users/ukaraoz/Work/Isogenie/IsogenieGenomes/define_guilds_forecosys_utils.R")

selected_traitsrulesgenes = select_traitsrulesgenes()
select_traits = selected_traitsrulesgenes[["select_traits"]]
select_rules = selected_traitsrulesgenes[["select_rules"]]
select_genes = selected_traitsrulesgenes[["select_genes"]]

microtrait.results.file = file.path(base, "IsogenieGenomes.microtraitresults.rds")
taxonomy.file = file.path(base, "genome_domain_phylogeny2.csv")
metadata.file = file.path(base, "sample2metadata.txt")
growthpredresults.file = file.path(base, "growthpredresults.rds")
genomelength.file = file.path(base, "IsogenieGenomes.length.txt")

#library(dplyr)
#load("growthpredresults.RData")
#growthpredresults = growthpredresults %>% tibble::rownames_to_column(var="id")%>%tibble::as_tibble()
#saveRDS(growthpredresults, file="growthpredresults.rds")

genomeset_results = readRDS(microtrait.results.file)
trait_matrixatgranularity3 = genomeset_results[["trait_matrixatgranularity3"]]
rule_matrix = genomeset_results[["rule_matrix"]]
hmm_matrix = genomeset_results[["hmm_matrix"]]
genomeset_growthpred = readRDS(growthpredresults.file) %>% dplyr::select(c(`id`, `mingentime`, `OGT`))
genomeset_metadata = build_metadata(taxonomy.file, metadata.file) %>% 
  dplyr::select(c(`genome`, `domain`, `phylum`, `habitat`, `depth`)) %>% 
  dplyr::rename(id = genome)
genomelength = read.table(genomelength.file, header = T, sep = "\t") %>% tibble::as_tibble()
# 
select_cols = c("id", "domain", "phylum", "habitat", "depth", "bp", "mingentime", "OGT")

ecosysguilds = trait_matrixatgranularity3 %>% 
  dplyr::left_join(rule_matrix, by = c("id" = "id")) %>% 
  dplyr::left_join(hmm_matrix, by = c("id" = "id")) %>% 
  dplyr::left_join(genomeset_growthpred, by = c("id" = "id")) %>%
  dplyr::left_join(genomeset_metadata, by = c("id" = "id"))  %>%
  dplyr::left_join(genomelength, by = c("id" = "id"))  %>%
  dplyr::select(c(select_cols,select_traits, select_rules, select_genes))
ecosysguilds = ecosysguilds %>%
  dplyr::mutate(ETC = 
                  case_when((get(selected_traitsrulesgenes[["ETC_traits"]][1])==1|
                    get(selected_traitsrulesgenes[["ETC_traits"]][2])==1|
                    get(selected_traitsrulesgenes[["ETC_traits"]][3])==1|
                    get(selected_traitsrulesgenes[["ETC_traits"]][4])==1)~as.integer(1), TRUE~as.integer(0))) %>%
  dplyr::mutate(heterotroph = 
                  case_when((get(selected_traitsrulesgenes[["heterotrophic_traits"]][1])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][2])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][3])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][4])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][5])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][6])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][7])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][8])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][9])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][10])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][11])==1|
                    get(selected_traitsrulesgenes[["heterotrophic_traits"]][12])==1)~as.integer(1), TRUE~as.integer(0))) %>%
  dplyr::mutate(fermentation = 
                  case_when((get(selected_traitsrulesgenes[["fermentation_traits"]][1])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][2])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][3])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][4])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][5])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][6])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][7])==1|
                    get(selected_traitsrulesgenes[["fermentation_traits"]][8])==1)~as.integer(1), TRUE~as.integer(0))) %>%
  dplyr::mutate(denitrification = 
                  case_when((get(selected_traitsrulesgenes[["denitrification_traits"]][1])==1|
                    get(selected_traitsrulesgenes[["denitrification_traits"]][2])==1|
                    get(selected_traitsrulesgenes[["denitrification_traits"]][3])==1|
                    get(selected_traitsrulesgenes[["denitrification_traits"]][4])==1)~as.integer(1), TRUE~as.integer(0)),
                ammonia_oxidation = 
                  case_when((get(selected_traitsrulesgenes[["NH4_NO2_oxidation_traits"]][1])==1)~as.integer(1), TRUE~as.integer(0)),
                nitrite_oxidation = 
                  case_when((get(selected_traitsrulesgenes[["NH4_NO2_oxidation_traits"]][2])==1)~as.integer(1), TRUE~as.integer(0)),
                nitrogenfixation = 
                  case_when((get(selected_traitsrulesgenes[["nitrogenfixation_traits"]][1])==1)~as.integer(1), TRUE~as.integer(0)),
                aerobic_diazotroph = 
                  case_when((get(selected_traitsrulesgenes[["nitrogenfixation_traits"]][1])==1&
                             ETC==1)~as.integer(1), TRUE~as.integer(0)),
                anaerobic_diazotroph = 
                  case_when((get(selected_traitsrulesgenes[["nitrogenfixation_traits"]][1])==1&
                             ETC==0)~as.integer(1), TRUE~as.integer(0)),
                methane_oxidation = 
                  case_when((get(selected_traitsrulesgenes[["methanemethanoloxidation_traits"]][1])==1&
                             ETC==1)~as.integer(1), TRUE~as.integer(0)),
                acetoclastic_methanogenesis = 
                  case_when((get(selected_traitsrulesgenes[["acetoclasticmethanogenesis_traits"]][1])==1&
                             (mcrA==1|mcrB==1|mcrG==1))~as.integer(1), TRUE~as.integer(0)),
                hydrogenotrophic_methanogenesis = 
                  case_when((get(selected_traitsrulesgenes[["hydrogenotrophicmethanogenesis_traits"]][1])==1&
                             (mcrA==1|mcrB==1|mcrG==1))~as.integer(1), TRUE~as.integer(0)),
                noETC = 
                  case_when(ETC==0~as.integer(1), TRUE~as.integer(0)),
                heterotroph_withETC = 
                  case_when((heterotroph==1&ETC==1)~as.integer(1), TRUE~as.integer(0))
                ) %>%
  dplyr::select(c(select_cols, `denitrification`, `ammonia_oxidation`, `nitrite_oxidation`, 
    `nitrogenfixation`, `aerobic_diazotroph`, `anaerobic_diazotroph`,
    `methane_oxidation`, `hydrogenotrophic_methanogenesis`, `acetoclastic_methanogenesis`,
    `heterotroph`, `heterotroph_withETC`, `ETC`, `noETC`, 
    select_rules, select_genes)) %>%
  dplyr::rename(`genome` = `id`)
write.table(ecosysguilds, file = file.path(base, "IsogenieGenomes.ecosysguilds.xls"), sep = "\t", 
            row.names = F, col.names = T, quote = F)

# some hmm level checks
# trait_check[["gene_check"]][["methanogenesis"]][["hydrogenotrophic"]]
# hmm_matrix_select = hmm_matrix %>% 
#   dplyr::mutate_if(is.factor, ~as.numeric(as.character(.x))) %>% 
#   dplyr::inner_join(genomeset_metadata, by = c("id" = "id")) %>% 
#   select(c(`id`, `domain`, `phylum`, trait_check[["gene_check"]][["methanogenesis"]][["hydrogenotrophic"]])) %>%
#   dplyr::mutate(sum = rowSums(across(where(is.numeric))))
# write.table(hmm_matrix_select, file = file.path(base, "IsogenieGenomes.hydrogenotrophic.xls"), sep = "\t", 
#             row.names = F, col.names = T, quote = F)





