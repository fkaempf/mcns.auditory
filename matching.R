library(arrow)
library(dplyr)
library(fafbseg)
library(malecns)
mcns_fw_matching = readRDS("/Users/fkampf/Documents/matching/10_2024_matching.RDS")
#baker_neurons2=fafbseg::add_celltype_info(baker_neurons, idcol='root_783')
#arrow::write_feather(baker_neurons2, sink = 'baker_neurons2.feather')

baker_neurons =read_feather("/Users/fkampf/Documents/matching/baker_neurons2.feather")
baker_neurons = baker_neurons %>% rename('neuron_name'='Neuron name')
baker_neurons <- baker_neurons%>% mutate(baker_type=sub("_[LR][0-9]*$", "", neuron_name))
keeps <- c("y", "a")
mba = mcns_body_annotations()

#stolen code from bella https://flyconnectome.slack.com/archives/D07QR4DKDQE/p1733247921466259
mcns_fw_f <- mcns_fw_matching %>% 
  select(query, match, conmatch) %>%
  mutate(query = gsub("_fw$", "", query)) %>%
  mutate(match = gsub("_mcns$", "", match)) %>% 
  mutate(conmatch = gsub("_mcns$", "", conmatch))

mcns_fw_f1 <- add_celltype_info(mcns_fw_f, idcol = "query", table = "info")

mcns_fw_f2 <- mcns_fw_f1 %>% 
  select(query, match, conmatch, cell_class, cell_type, ito_lee_hemilineage) %>% 
  rename("fw_hl" = "ito_lee_hemilineage", "fw_type" = "cell_type", "fw_class" = "cell_class")

mcns_fw_f2$match <- as.integer(mcns_fw_f2$match)
mcns_fw_f2$conmatch <- as.integer(mcns_fw_f2$conmatch)

# First, merge to get the `match_type`
mcns_fw_f3 <- mcns_fw_f2 %>%
  left_join(mba %>% select(bodyid, type,itolee_hl), by = c("match" = "bodyid")) %>%
  rename(match_type = type)%>% # Rename the column for match_type
  rename(match_itolee_hl = itolee_hl)  # Rename the column for hemilineage

# Second, merge to get the `conn_match_type`
mcns_fw_f3 <- mcns_fw_f3 %>%
  left_join(mba %>% select(bodyid, type,itolee_hl), by = c("conmatch" = "bodyid")) %>%
  rename(conn_match_type = type) %>% # Rename the column for conn_match_type
  rename(conn_match_itolee_hl = itolee_hl)   # Rename the column for conn_match_type hemilineage


#get unique fw_types
mcns_fw_f4<- mcns_fw_f3 %>% 
  group_by(fw_type) %>% 
  summarise(
    conn_match_type_combined = ifelse(
      length(unique(na.omit(conn_match_type))) == 0, NA, 
      paste(unique(na.omit(conn_match_type)), collapse = ", ")
    ),
    match_type_combined = ifelse(
      length(unique(na.omit(match_type))) == 0, NA, 
      paste(unique(na.omit(match_type)), collapse = ", ")
    )
  )

mcns_fw_f3 %>% 
  count(fw_type, conn_match_type, match_type) %>% 
  group_by(fw_type) %>%
  summarise(match_type=paste(na.omit(match_type), collapse = ", "),
            conn_match_type=paste(na.omit(conn_match_type), collapse = ", "))
  


#get unique cell types that need to be matched from the baker neurons
unique_cell_types = unique(baker_neurons$cell_type)
baker2mcns.matches.df <- tibble(
  cell_type_fw = unique_cell_types) %>%
  mutate(type = NA)
print(nrow(baker2mcns.matches.df))

#look for matches that already exist in the male body annotation (mba)


# option 1 ...
mba.match_exist <-mba %>% 
  filter((flywire_type %in% unique_cell_types & !is.na(flywire_type)) |
           (!is.na(type) & type %in% unique_cell_types) )%>%
  filter(!is.na(flywire_type))


# ... option2
# make a new copy of mba (malecns annotations) with one column called "ptype"
# which has the flywire type if that exists or the "type" column if not
mba.pf=mba
mba.pf$ptype=mcns_predict_type(mba, prefer.foreign = T)
mba.match_exist2 <-mba.pf %>% 
  filter(ptype %in% baker2mcns.matches.df$cell_type_fw & !is.na(ptype))
nrow(mba.match_exist3)

# select malecns neurons that we can match to Baker neurons
mba.match4 <- mba.match_exist2 %>% 
  select(bodyid, flywire_type) %>% 
  left_join(baker_neurons %>% distinct(cell_type, .keep_all = T), c("flywire_type"="cell_type"))
#  mutate(baker_type=sub("_[LR][0-9]*$", "", neuron_name)) %>% 
#  mutate(baker_type=sub("_[0-9]+$", "", baker_type))

expected_secondary <- c("A1", "A2", "AVLP_pr23", "WV-WV", "GF", "B1", "B2")

# let's look at second order neurons (that receive direct sensory input)
mba.match5<- mba.match4 %>% 
  filter(grepl(paste(expected_secondary, collapse = "|") , neuron_name))%>%
  mutate(baker_type=sub("_[LR][0-9]*$", "", neuron_name)) %>% 
  distinct(bodyid,.keep_all = T) %>%
  arrange(baker_type)
#85













