library(arrow)
library(dplyr)
library(fafbseg)
library(malecns)
mcns_fw_matching = readRDS("/Users/fkampf/Documents/matching/10_2024_matching.RDS")
baker_neurons =read_feather("/Users/fkampf/Documents/matching/baker_neurons2.feather")
keeps <- c("y", "a")
mba = mcns_body_annotations()

#stolen code
mcns_fw_f <- mcns_fw_matching %>% select(query, match, conmatch) %>%
  mutate(query = gsub("_fw$", "", query)) %>% mutate(match = gsub("_mcns$", "", match)) %>% mutate(conmatch = gsub("_mcns$", "", conmatch))

mcns_fw_f1 <- add_celltype_info(mcns_fw_f, idcol = "query", table = "info")

mcns_fw_f2 <- mcns_fw_f1 %>% select(query, match, conmatch, cell_class, cell_type, ito_lee_hemilineage) %>% rename("fw_hl" = "ito_lee_hemilineage", "fw_type" = "cell_type", "fw_class" = "cell_class")

mcns_fw_f2$match <- as.integer(mcns_fw_f2$match)
mcns_fw_f2$conmatch <- as.integer(mcns_fw_f2$conmatch)


# First, merge to get the `match_type`
mcns_fw_f3 <- mcns_fw_f2 %>%
  left_join(mba %>% select(bodyid, type), by = c("match" = "bodyid")) %>%
  rename(match_type = type)  # Rename the column for match_type

mcns_fw_f3 <- mcns_fw_f3 %>%
  left_join(mba %>% select(bodyid, itolee_hl), by = c("match" = "bodyid")) %>%
  rename(match_itolee_hl = itolee_hl)  # Rename the column for match_type

# Second, merge to get the `conn_match_type`
mcns_fw_f3 <- mcns_fw_f3 %>%
  left_join(mba %>% select(bodyid, type), by = c("conmatch" = "bodyid")) %>%
  rename(conn_match_type = type)  # Rename the column for conn_match_type

mcns_fw_f3 <- mcns_fw_f3 %>%
  left_join(mba %>% select(bodyid, itolee_hl), by = c("conmatch" = "bodyid")) %>%
  rename(conn_match_itolee_hl = itolee_hl)  # Rename the column for conn_match_type

unique_cell_types = unique(baker_neurons$cell_type)


