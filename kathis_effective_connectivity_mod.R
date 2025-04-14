library(neuprintr)
library(arrow)
library(malecns)
library(bit64)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(cowplot)

#get the whole maleCNS matrix, these first steps will take quite some time
mdf = mcns_dvid_annotations()
query.ids = as.integer64(mdf$bodyid[!(mdf$status %in% c("Orphan","PRT Orphan", "RT Orphan",
                                                        "Unimportant","Orphan hotknife","Out of scope",
                                                        "Glia","Orphan-artifact",NA))])
query_info = mcns_neuprint_meta(query.ids)
query_info = query_info[!(query_info$status %in% c("Orphan",NA)),]
query.ids2 = na.omit(as.character(query_info$bodyid))
#clean up along the way
rm(query.ids,query_info)

#get the flat connectome
conn_feather <-  arrow::open_dataset('/Users/fkampf/Downloads/connectome-weights-2023-11-15-f20f31-unlocked-minconf-0.5-primary-only.feather', format = 'feather')
conn_feather.inmem <- conn_feather %>%
  select(body_pre, body_post, weight) %>%
  collect()

# NB make sure that we have a square matrix with all ids (even if some have no inputs/outputs)
all_ids=union(unique(conn_feather.inmem$body_pre), unique(conn_feather.inmem$body_post))
np_adj <- coconat::partner_summary2adjacency_matrix(conn_feather.inmem, inputcol = 'body_pre', outputcol = 'body_post', sparse = T,
                                                    inputids = all_ids, outputids = all_ids,
                                                    standardise_input=F)

#' Efficient column and row scaling for sparse matrices
#'
#' @details Conventionally in connectivity matrices, rows are upstream input
#'   neurons and columns are downstream output neurons. \code{colScale} will
#'   therefore normalise by the total input onto each downstream neuron while
#'   \code{colScale} will normalise by the total output from each upstream
#'   neuron.
#'
#'   These functions will work with dense inputs but are less efficient than
#'   their base R cousins (e.g. \code{scale}). As an example, normalising a 4x4k
#'   subset of VNC connectivity data took \itemize{
#'
#'   \item colScale sparse 3.5ms
#'
#'   \item colScale dense 1460ms
#'
#'   \item scale dense 365ms
#'
#'   }
#'
#'   In conclusion sparse matrices with colScale/rowScale are the way to go!
#' @param A a sparse (or dense) matrix
#' @param na.rm
#'
#' @return A scaled sparse matrix
#' @export
#'
#' @seealso
#' \url{https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r}
#'
#' @examples
#' # box::use(mancfuns/manc)
#' library(Matrix)
#' set.seed(42); A <- Matrix(rbinom(100,10,0.05), nrow = 10)
#' manc$colScale(A)
#' manc$rowScale(A)
#' geomscalea=sqrt(manc$colScale(A)*manc$rowScale(A))
colScale <- function(A, na.rm=TRUE) {
  scalefac = 1 / Matrix::colSums(A)
  if(na.rm) scalefac[!is.finite(scalefac)]=0
  B = A %*% Matrix::Diagonal(x = scalefac)
  B
}

#now calculate the input percent using this funtion
np_adj_per_in <- colScale(np_adj)
colnames(np_adj_per_in) = colnames(np_adj)

#restrict matrix to real neuron bodies
np_adj <- np_adj[rownames(np_adj) %in% query.ids2,colnames(np_adj) %in% query.ids2]
np_adj_per_in <- np_adj_per_in[rownames(np_adj_per_in) %in% query.ids2,colnames(np_adj_per_in) %in% query.ids2]

#clean up
rm(query.ids2,conn_feather.inmem,conn_feather,mdf)

#set the threshold and number of layers
desired_layers = 3

#get your input IDs and info from clio
all_info = mcns_body_annotations()

cell_types <- c('CB1076','CB1078', 'CB3710', 'CB2521','CB1542', 'CB1038', 'CB1427', 'CB2556', 'CB2380')
mba<-mcns_body_annotations()
mba.type.coalesced <-  mba %>%
  mutate(type = coalesce(type, 
                         flywire_type, 
                         hemibrain_type, 
                         manc_type))
input_ids <-mba.type.coalesced %>% 
  filter(flywire_type %in% cell_types)%>%
  pull(bodyid)%>%
  unique()

input_info = mba[mba$bodyid %in% input_ids,]
input_info_sel <- input_info[,c("bodyid","class","group","type","instance","soma_side", "root_side", "superclass")]
input_info_sel$soma_or_root_side = ifelse(!is.na(input_info_sel$soma_side)&input_info_sel$soma_side!="",input_info_sel$soma_side,input_info_sel$root_side)


output_ids <- as.character(c(57956,30261,27091,18554,17337,13495,17918,19136,515347,15035,16719,
                             17089,17589,31586,12442,18276,17098,20597,267996,24537,21349,79094,32788,27072,
                             16559,43592,21586,13655,22744,20803,15324,84181,518589,18830,19894,556999,30669,
                             19818,513129,517178,559181,35694,27507,53355,532582,17867,555488,514909,24072,
                             17618,17851,13790,22346,17434,10930,70011,19948,29992,513515,513651,517653,
                             516954,20392,10765,35993,539247,21507,37651,30273,16805,20117,555751,13562,
                             21033,29043,20068,87497,522357,27330,12889,75407,34082,27368,41366,26543,19600,
                             24874,33318,529875,41518,24060,23968,23558,518253,16541,528986,555811,31905,
                             514077,55707,33009,17841,19109,33503,517337,519518,25542,22093,521453,521449,
                             531383,15876,46660,518038,30391,17702,515990,37522,18861,33435,47325,29462,
                             14492,557035,513652,17487,18017,55136,30558,19200,19236,524297,17470,29483,
                             30455,18990,50692,37977,29143,26989,28094,33519,522419,71385,26746,23319,22270))

output_info = all_info[all_info$bodyid %in% output_ids,]
output_info_sel <- output_info[,c("bodyid","class","group","type","instance","soma_side", "root_side", "superclass")]
output_info_sel$soma_or_root_side = ifelse(!is.na(output_info_sel$soma_side)&output_info_sel$soma_side!="",output_info_sel$soma_side,output_info_sel$root_side)

#precalculate all layers for all input neurons, this step will take a while depending on how many input neurons you are using
#set up a list of vectors for each iteration through the network and make the first vector considering only your starting neurons
layer_vec = paste0("v",seq(from=0,to=desired_layers))
layer_list <- vector("list", length(layer_vec))
names(layer_list) <- layer_vec
layer_list[["v1"]] = np_adj_per_in[rownames(np_adj_per_in) %in% input_info_sel$bodyid,]
for (layer in 2:desired_layers) {
  layer_list[[layer+1]] <-  layer_list[[layer]] %*% np_adj_per_in
}

#now pull out data from multiplied matrices for input neurons, this can be based on a lot of things, choosing cell_sub_class as example
input_type <- unique(input_info_sel$type)
for (s in 1:length(input_type)) {
  assign(input_type[s], na.omit(unique(input_info_sel$bodyid[grepl(input_type[s],input_info$type)])))
}
#for output neurons, this can be based on a lot of things, choosing type as example
output_order <- sort(unique(output_info$type))


in_out_scores <- setNames(data.frame(matrix(ncol = 16, nrow = 0)),c("bodyid","conn_score","ds_layer","class","group","type","instance","soma_side","rootSide",
                                                                    "superclass","soma_or_root_side","input","in_side","ipsi_contra")) #, "clone"
for (n in 1:length(input_type)) {
  startn_r = as.character(na.omit(input_info_sel$bodyid[input_info_sel$type == input_type[n] & input_info_sel$soma_or_root_side == "R"]))
  startn_l = as.character(na.omit(input_info_sel$bodyid[input_info_sel$type == input_type[n] & input_info_sel$soma_or_root_side == "L"]))
  if (length(startn_r)>0 & length(startn_l)>0) {
    input = input_type[n]
    
    #set up a list of vectors for each iteration through the network
    norm_list_l <- vector("list", length(layer_vec))
    names(norm_list_l) <- layer_vec
    norm_list_r <- vector("list", length(layer_vec))
    names(norm_list_r) <- layer_vec
    #get data from precalculated matrices and normalize
    for (layer in 1:desired_layers) {
      if (length(startn_l)>1) {
        temp_layer_conn_str = Matrix::colSums(layer_list[[layer+1]][na.omit(match(startn_l, rownames(layer_list[[layer+1]]))),])
      } else {
        temp_layer_conn_str = layer_list[[layer+1]][na.omit(match(startn_l, rownames(layer_list[[layer+1]]))),]
      }
      norm_list_l[[layer+1]] <- setNames(as.vector(temp_layer_conn_str/mean(temp_layer_conn_str[temp_layer_conn_str>0])),
                                         colnames(np_adj_per_in))
      if (length(startn_r)>1) {
        temp_layer_conn_str = Matrix::colSums(layer_list[[layer+1]][na.omit(match(startn_r, rownames(layer_list[[layer+1]]))),])
      } else {
        temp_layer_conn_str = layer_list[[layer+1]][na.omit(match(startn_r, rownames(layer_list[[layer+1]]))),]
      }
      norm_list_r[[layer+1]] <- setNames(as.vector(temp_layer_conn_str/mean(temp_layer_conn_str[temp_layer_conn_str>0])),
                                         colnames(np_adj_per_in))
    }
    
    #now get the highest connectivity score for the input to each output type
    out_r_scores_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("bodyid","conn_score","ds_layer"))
    for (layer in 2:desired_layers) {
      out_r_scores <- norm_list_r[[layer]][names(norm_list_r[[layer]])%in% output_ids]
      df = data.frame(bodyid=as.double(names(out_r_scores)), conn_score=out_r_scores, ds_layer = layer-1, row.names=NULL)
      out_r_scores_df = rbind(out_r_scores_df, df)
    }
    out_r_scores_df %>%
      group_by(bodyid) %>%
      dplyr::slice(which.max(conn_score)) -> out_r_scores_df_max
    
    out_r_scores_df_max_info <- left_join(out_r_scores_df_max, output_info_sel, by = "bodyid")
    out_r_scores_df_max_info$input <- input
    out_r_scores_df_max_info$in_side <- "R"
    
    out_l_scores_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("bodyid","conn_score","ds_layer"))
    for (layer in 2:desired_layers) {
      out_l_scores <- norm_list_l[[layer]][names(norm_list_l[[layer]])%in% output_ids]
      df = data.frame(bodyid=as.double(names(out_l_scores)), conn_score=out_l_scores, ds_layer = layer-1, row.names=NULL)
      out_l_scores_df = rbind(out_l_scores_df, df)
    }
    out_l_scores_df %>%
      group_by(bodyid) %>%
      dplyr::slice(which.max(conn_score)) -> out_l_scores_df_max
    
    out_l_scores_df_max_info <- left_join(out_l_scores_df_max, output_info_sel, by = "bodyid")
    out_l_scores_df_max_info$input <- input
    out_l_scores_df_max_info$in_side <- "L"
    
    out_scores <- rbind(out_r_scores_df_max_info,out_l_scores_df_max_info)
    out_scores$ipsi_contra <- NA
    out_scores$ipsi_contra[out_scores$soma_or_root_side == out_scores$in_side] <- "ipsi"
    out_scores$ipsi_contra[out_scores$soma_or_root_side != out_scores$in_side] <- "contra"
    
    in_out_scores <- rbind(in_out_scores,out_scores)
  }
}

in_out_scores %>% group_by(input,type,ipsi_contra) %>%
  summarise_at(vars("conn_score", "ds_layer"), mean) -> in_out_scores_mean
#cluster result on ipsi connectivity
in_out_scores_mean %>% filter(ipsi_contra=="ipsi") -> in_out_scores_mean_ipsi
in_out_scores_mean_ipsi_m <- tidyr::pivot_wider(in_out_scores_mean_ipsi[,c("input","conn_score","type")], names_from = "input", values_from = "conn_score")
in_out_scores_mean_ipsi_m <- as.matrix(in_out_scores_mean_ipsi_m[, -1]) # -1 to omit categories from matrix

clust <- hclust(dist(t(in_out_scores_mean_ipsi_m)))
my_cols_fun <- colorRampPalette(c("#3E1F4B","#8B008B","#F2F2F2"))

#plot as heatmap for all input to output neurons
figure_heatmap <- ggplot(as_tibble(in_out_scores_mean), aes(x = factor(type, level = output_order), y = input)) +
  geom_point(aes(col = ds_layer, size = conn_score), shape = 15) +
  theme_minimal() +
  theme(
    legend.position = 'right',
    text = element_text(color = 'grey40')
  ) +
  scale_size_area(max_size = 5) +
  scale_colour_gradientn(colours = my_cols_fun(desired_layers), limits=c(1, desired_layers)) +
  scale_y_discrete(limits = colnames(in_out_scores_mean_ipsi_m)[clust$order]) +
  guides(colour = guide_legend(override.aes = list(size=10)))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 14, family = "Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) #+
#facet_grid(. ~ factor(ipsi_contra, levels=c("ipsi","contra")), scales = "fixed")

#boxplots for the same data above
for (nth in 1:length(input_types)) {
  in_out_scores_type <- in_out_scores[in_out_scores$input == input_types[nth] & in_out_scores$type %in% output_order,]
  
  type_score <- ggplot(in_out_scores_type, aes(x=factor(type, level = output_order), y=conn_score)) +
    geom_boxplot(fill="#8B008B") + ylab("Score") + xlab("output neurons") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_line(colour = "grey", linewidth = 0.2)) + ylim(0, ceiling(max(in_out_scores_type$conn_score))) +
    facet_grid(. ~ factor(ipsi_contra, levels=c("ipsi","contra")), scales = "fixed")
  
  type_layer <- ggplot(in_out_scores_type, aes(x=factor(type, level = output_order), y=ds_layer)) +
    geom_boxplot(fill="#8B008B") + ylab("Layer") + xlab("output neurons") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_line(colour = "grey", linewidth = 0.2)) + ylim(1, desired_layers) +
    facet_grid(. ~ factor(ipsi_contra, levels=c("ipsi","contra")), scales = "fixed")
  
  legend <- get_legend(
    type_score +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  title <- ggdraw()+
    draw_label(paste0(input_types[nth]," connectivity to leg MNs"), fontface='bold')
  
  figure_nwh <- plot_grid(type_score + theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position="none",plot.margin = margin(0.1,0.1,0,0.3, "cm")),
                          type_layer + theme(legend.position="none",strip.background = element_blank(),strip.text.x = element_blank(),plot.margin = margin(0.1,0.1,0,0.3, "cm")),
                          ncol = 1, common.legend = FALSE,
                          rel_heights = c(1.3, 0.8), align = "v", vjust=0, axis = "rlbt")
  figure2_nwh <- plot_grid(title, legend, nrow=1, rel_widths =c(1,1))
  figure3_nwh <- plot_grid(figure2_nwh, figure_nwh, ncol=1, rel_heights =c(0.1,1))
  ggsave(paste0("output/",input_types[nth], "_combined_plot.pdf"), figure3_nwh, dev=cairo_pdf, width=4715, height=3295, units= "px")
}


