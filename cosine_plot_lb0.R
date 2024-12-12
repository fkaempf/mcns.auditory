library(fafbseg)
library(coconatfly)
library(arrow)
library(dplyr)

baker_neurons =read_feather("/Users/fkampf/Documents/matching/baker_neurons2.feather")
baker_neurons = baker_neurons %>% rename('neuron_name'='Neuron name')
lb0neurons=cf_meta(cf_ids(malecns='/itoleeHl:LB0.*', flywire = baker_neurons %>% filter(grepl("WV", neuron_name)) %>% pull(root_783) ))
lb0neurons=cf_meta(cf_ids(malecns='/itoleeHl:LB0.*', flywire = baker_neurons$root_783 ))

lb0neurons2=left_join(lb0neurons, baker_neurons %>% select(root_783, neuron_name), by=c("id"="root_783"))
lb0neurons3 <- lb0neurons2 %>% 
  mutate(instance=case_when(
    is.na(neuron_name) ~ instance,
    T ~ neuron_name
  ))

lb0neurons3 %>% 
  with(cf_cosine_plot(key, labRow = instance, interactive = T))
