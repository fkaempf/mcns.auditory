library(tidyr)
library(dplyr)
library(malecns)
library(arrow)
library(fafbseg)
library(coconatfly)

baker_neurons = read_feather("/Users/fkampf/Documents/matching/baker_neurons2.feather")
baker_neurons = baker_neurons %>% rename('neuron_name'='Neuron name')

mba <- mcns_body_annotations()
lb0_male <- mba%>%filter(itolee_hl == "LB0") %>% select(bodyid) %>% pull() %>% as.character()
baker_fw <- baker_neurons %>% filter(grepl("WV", neuron_name)) %>% select(root_783) %>% pull()

lb0neurons_fw_meta=cf_meta(cf_ids(flywire = baker_fw))
lb0neurons_fw_meta <- lb0neurons_fw_meta %>% select(id, type, ito_lee_hemilineage, malecns_type, instance) %>% rename(bodyid = "id", itolee_hl = "ito_lee_hemilineage", flywire_type = "type")

lb0neurons_mcns_meta <- mba %>% filter(itolee_hl == "LB0") %>% select(bodyid, type, itolee_hl, flywire_type, instance) %>% rename(malecns_type = "type")

lb0neurons_all <- rbind(lb0neurons_fw_meta, lb0neurons_mcns_meta)

lb0neurons_all_info <- lb0neurons_all %>%
  left_join(baker_neurons, by =c("bodyid" = "root_783")) %>%
  mutate(instance = case_when(
    is.na(neuron_name) ~ instance,
    T ~ neuron_name
  ))

lb0neurons_all_info <- lb0neurons_all_info %>% filter(!is.na(flywire_type)) %>%
  mutate(flywire_type = ifelse(is.na(flywire_type), malecns_type, flywire_type))

lb0neurons_all_info$neuron_name <- sub("_.*", "", lb0neurons_all_info$neuron_name)

flywire_lookup <- lb0neurons_all_info %>%
  filter(!is.na(neuron_name)) %>% # Use only rows with valid neuron_name
  distinct(flywire_type, neuron_name) # Remove duplicates

flywire_lookup <- flywire_lookup %>%
  group_by(flywire_type) %>%
  summarise(neuron_name = first(neuron_name))

lb0neurons_all_info <- lb0neurons_all_info %>%
  left_join(flywire_lookup, by = "flywire_type", suffix = c("", "_lookup")) %>%
  mutate(neuron_name = coalesce(neuron_name, neuron_name_lookup)) %>%
  select(-neuron_name_lookup)

lb0neurons_all_info1 <- lb0neurons_all_info %>%
  filter(!is.na(neuron_name)) %>%
  select(bodyid, flywire_type, malecns_type, neuron_name, itolee_hl)

googlesheets4::write_sheet(data = lb0neurons_all_info1, ss= "https://docs.google.com/spreadsheets/d/1PMcHcy99lgJi_2Sk4msTwKfTs7NcK7fVGz40qLiXd1k/edit?usp=sharing")