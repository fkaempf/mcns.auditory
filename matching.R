library(arrow)
library(dplyr)
library(fafbseg)
library(malecns)
mcns_fw_matching = readRDS("/Users/fkampf/Documents/matching/10_2024_matching.RDS")
#baker_neurons2=fafbseg::add_celltype_info(baker_neurons, idcol='root_783')
#arrow::write_feather(baker_neurons2, sink = 'baker_neurons2.feather')

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
  left_join(mba %>% select(bodyid, type,itolee_hl), by = c("match" = "bodyid")) %>%
  rename(match_type = type)%>% # Rename the column for match_type
  rename(match_itolee_hl = itolee_hl)  # Rename the column for hemilineage

# Second, merge to get the `conn_match_type`
mcns_fw_f3 <- mcns_fw_f3 %>%
  left_join(mba %>% select(bodyid, type,itolee_hl), by = c("conmatch" = "bodyid")) %>%
  rename(conn_match_type = type) %>% # Rename the column for conn_match_type
  rename(conn_match_itolee_hl = itolee_hl)   # Rename the column for conn_match_type hemilineage

unique_cell_types = unique(baker_neurons$cell_type)
baker2mcns_map  <- data.frame(
  cell_type = unique_cell_types
)

#first take where flywire and mcns type already match
types.match.TRUE.bool <- replace(mba$type==mba$flywire_type,is.na(mba$type==mba$flywire_type),FALSE)
types.match.TRUE.df <- mba %>% filter(types.match.TRUE.bool)
types.match.TRUE.names <- mba %>% filter(types.match.TRUE.bool) 
types.match.FALSE.names <- mba %>% filter(!types.match.TRUE.bool) %>% select(type) %>%unique()



mcns_fw_f3.both <- mcns_fw_f3 %>% filter(!is.na(match_type)) %>% filter(match_type == conn_match_type)
mcns_fw_f3.nblast <- mcns_fw_f3 %>% filter(!is.na(match_type))
mcns_fw_f3.conn <- mcns_fw_f3 %>% filter(is.na(match_type)&!is.na(conn_match_type))

mcns_fw_f3.exist.matches <- types.match.TRUE.df %>% select(flywire_type,type)

before_matching_n <- nrow(types.match.FALSE.names)
mcns_fw_f3.both.matches <- mcns_fw_f3.both %>% filter(match_type %in% types.match.FALSE.names$type) %>% select(fw_type,match_type) %>% filter(!is.na(fw_type)) %>% unique()
types.match.FALSE.names <- types.match.FALSE.names %>% filter(!(type %in% mcns_fw_f3.both.matches$fw_type))

mcns_fw_f3.nblast.matches <- mcns_fw_f3.nblast %>% filter(match_type %in% types.match.FALSE.names$type) %>% select(fw_type,match_type) %>% filter(!is.na(fw_type)) %>% unique()
types.match.FALSE.names <- types.match.FALSE.names %>% filter((type %in% mcns_fw_f3.nblast.matches$fw_type))

mcns_fw_f3.conn.matches <- mcns_fw_f3.conn %>% filter(conn_match_type %in% types.match.FALSE.names$type) %>% select(fw_type,conn_match_type) %>% filter(!is.na(fw_type)) %>% unique()
types.match.FALSE.names <- types.match.FALSE.names %>% filter(!(type %in% mcns_fw_f3.conn.matches$fw_type))

print('Number of neurons matched')
print((before_matching_n - nrow(types.match.FALSE.names))/before_matching_n)

baker_neurons_matched_exist = baker_neurons %>% left_join(mcns_fw_f3.exist.matches%>%rename('match_type_exist'='type'),by=c('cell_type'='flywire_type')) %>% filter(is.na(match_type))
baker_neurons_matched_both = baker_neurons %>% left_join(mcns_fw_f3.both.matches%>%rename('match_type_both'='match_type'),by=c('cell_type'='fw_type')) %>% filter(is.na(match_type))
baker_neurons_matched_nblast = baker_neurons %>% left_join(mcns_fw_f3.nblast.matches%>%rename('match_type_nblast'='match_type'),by=c('cell_type'='fw_type')) %>% filter(is.na(match_type))
baker_neurons_matched_conn = baker_neurons_matched %>% left_join(mcns_fw_f3.conn.matches%>%rename('match_type_conn'='match_type'),by=c('cell_type'='fw_type')) %>% filter(is.na(match_type))


#what I thought is correct
test <- baker2mcns_map %>% 
  left_join(mcns_fw_f3%>%
              select(fw_type,match_type,conn_match_type),by=c('cell_type'='fw_type')) %>%
  unique() %>%
  filter(!is.na(cell_type)) %>%
  filter(!is.na(match_type)&!is.na(conn_match_type)) %>%
  group_by(cell_type)

baker_neurons_matched <-baker_neurons %>% left_join(test,by = c('cell_type'='cell_type'))