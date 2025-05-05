fafbseg::flywire_connectome_data("syn")
library(malecns)
library(coconatfly)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(stringr)
library(igraph)
library(tidygraph)
library(ggraph)
synapse_threshold=5

prepare_data <- function(dataset.name,partner.df,post_or_pre){
  if (grepl(post_or_pre, 'post')) {post_or_pre='post_key'}
  else {post_or_pre='pre_key'}
  
  data <- cf_meta(partner.df%>%
                    filter(dataset==dataset.name)%>%
                    pull(post_or_pre)%>% #should be post right now
                    unique())%>%
    mutate(id=as.character(id))%>%
    rename_with(~ str_replace_all(.x, c(
      "flywireType" = "flywire",
      "flywire_type" = "flywire",
      "hemibrainType" = "hemibrain",
      "hemibrain_type" = "hemibrain",
      'malecns_type' = 'malecns',
      'malecns_type' = 'malecns',
      "type" = 'type2',
      'bodyid'='id'
    )))%>%
    select(any_of(c('id',
                    'type1',
                    'type2',
                    'flywire',
                    'hemibrain',
                    'malecns')))
  
  return(data)
}

make_edges <- function(data,sink_no=1){
  data %>% 
    pivot_longer(!id,
                 names_to = 'dataset',
                 values_to = "type")%>%
    select(id,type)%>%
    filter(!is.na(type))%>%
    mutate(sink=sink_no)}







match_cell_types <-function(partner.df,post_or_pre='post',flywire2mcns=T){
  #format data
  if (flywire2mcns){  
    data1 <- prepare_data('flywire', partner.df,post_or_pre)
    data2 <- prepare_data('malecns', partner.df,post_or_pre)}
  else {  
    data1 <- prepare_data('malecns', partner.df,post_or_pre)
    data2 <- prepare_data('flywire', partner.df,post_or_pre)}
  
  

  
  #concat data
  data.concat <- rbind(data1%>%
                         rename('type'='type2')%>%
                         select(id,type),
                       data2%>%
                         rename('type'='type2')%>%
                         select(id,type))
  
  #make edges
  edges1 <- make_edges(data1,'1')
  edges2 <- make_edges(data2,'2')
  edges<-rbind(edges1,edges2)
  
  #make nodes from edges
  nodes<-  edges %>% 
    pivot_longer(cols = c('id','type'), values_to = "name",names_to = 'type')%>%
    mutate(node_type=if_else(type=='type',
                             'cell_type',
                             str_c("sink_", sink)))%>%
    select(name, node_type)%>%
    distinct(name,.keep_all = T)
  
  #create graph
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  
  # Process graph distances
  graph.distances <- data.frame(distances(graph))
  colnames(graph.distances)<- nodes$name
  rownames(graph.distances)<- nodes$name
  
  graph.distances.sink2sink <- graph.distances[nodes$node_type=='sink_1',nodes$node_type=='sink_2']
  rownames(graph.distances.sink2sink)<- nodes$name[nodes$node_type=='sink_1']
  colnames(graph.distances.sink2sink)<- nodes$name[nodes$node_type=='sink_2']
  
  column.indices.shortest.path <- sapply(rownames(graph.distances.sink2sink), function(x) {
    temp <- graph.distances.sink2sink[x, , drop = FALSE]
    
    if(length(unique(temp)) == 1) {
      return(NA) 
    } else {
      min_col <- which(temp == min(temp), arr.ind = TRUE)[1,][2]
      
      if (temp[1, min_col] == Inf) {
        return(NA)
      }
      
      return(min_col)
    }
  })
  column.indices.shortest.path <- setNames(column.indices.shortest.path, rownames(graph.distances.sink2sink))
  
  matched.type.list<-list()
  
  
  for (x in rownames(graph.distances.sink2sink)) {
    
    
    
    if (is.na(column.indices.shortest.path[x])) {
      
      print('No matched type')
      matched.type = data.concat%>%filter(id==x)%>%pull(type)
      
    }
    else{
      node1 <- x
      node2 <- colnames(graph.distances.sink2sink)[column.indices.shortest.path[x]]
      
      # Get the shortest path once
      path <- shortest_paths(graph, from = node1, to = node2)
      
      matched.type <- V(graph)$name[path$vpath[[1]][length(path$vpath[[1]]) - 1]]
      if(length(matched.type)==0){readline()
        
      }
      writeLines(paste(paste("\nFlywire_id", node1), 
                       paste("\nmalecns_id", node2),
                       paste("\nmatched_type", matched.type),
                       '\n'))
    }
    
    matched.type.list <- append(matched.type.list, list(matched.type))
    
  }
  matched.join.df <-  data.frame(id = rownames(graph.distances.sink2sink), match.type = as.character(matched.type.list))
  
  result<- rbind(data1%>%
                   left_join(matched.join.df,by='id')%>%
                   select(id,match.type),
                 data2%>%rename('match.type'='type2')%>%
                   select(id,match.type))
  return(result)
  
}


cell_types <- c('CB1076', 'CB3710', 'CB2521', 'CB1425', 'CB1078', 'CB1038', 'CB1427', 'CB2556','CB2380','CB2556')  #ALL
cf_list <- lapply(cell_types, function(ct) {
  cf_ids(malecns = paste0("flywireType:", ct), flywire = ct)
})
all.auditory.output <- cf_partners(do.call(c, cf_list), partners = 'out', threshold = synapse_threshold)%>%
  mutate(post_id=as.character(post_id))

result.f2m <- match_cell_types(all.auditory.output, 'post',flywire2mcns=T)%>%rename('match.type.f2m'='match.type')
result.m2f <- match_cell_types(all.auditory.output, 'post',flywire2mcns=F)%>%rename('match.type.m2f'='match.type')


all.auditory.output.match <- all.auditory.output %>%left_join(result.f2m,by=c('post_id'='id'))

