#' Create network visualization
#' 
#' @importFrom magrittr "%>%"
#' @param node_df data.frame of nodes in network
#' @param edge_df data.frame of edges in network
#' @return An igraph object (see package \code{igraph} for more details)
#' @details 
#' `node_df` must have a column called `id` that contains the node identifiers. `edge_df` must have at least
#' two columns labeled `from` and `to` respectively. Optional columns in `node_df` are `label`, `size`,
#' `color`, and `shape`
#' @examples
#' nodes <- data.frame(id=LETTERS)
#' edges <- data.frame(from=LETTERS[sample(1:26, 10)], to=LETTERS[sample(1:26, 30, replace=T)])
#' make_a_network(nodes, edges)
#' @export make_a_network

make_a_network <- function(node_df, edge_df,
                           directed=T, layout=layout_nicely,
                           plot=T, 
                           detect_communities=T,
                           detect_cliques=F){
  
  require(igraph)
  
  graph <- igraph::graph_from_data_frame(d=edge_df, vertices=node_df,
                                directed=directed)
  
  if(is.null(node_df[["label"]])){
    vertex_label <- node_df[["id"]] %>% as.character()
  } else {
    vertex_label <- node_df[["label"]] %>% as.character()
  }
  
  if(is.null(node_df[["size"]])){
    vertex_size <- 7
  } else {
    vertex_size <- node_df[["size"]] %>% 
                        as.numeric() %>% 
                        abs() %>% 
                        magrittr::divide_by(max(.)) %>%
                        magrittr::multiply_by(10)
  }
  
  if(is.null(node_df[["color"]])){
    vertex_color <- "#67a9cf"
  } else {
    vertex_color <- node_df[["color"]]
  }
  
  if(is.null(node_df[["shape"]])){
    vertex_shape <- "circle"
  } else {
    vertex_shape <- node_df[["shape"]]
  }
  
  if(detect_communities){
    ceb <- igraph::cluster_edge_betweenness(graph)
    vertex_color <- RColorBrewer::brewer.pal(12, "Set3") %>%
                      magrittr::extract((ceb$membership %% 12) + 1)
  }
  
  if(plot){
    graph_plot <- igraph::plot.igraph(graph, 
                                      vertex.label.dist=1,
                                      vertex.label=vertex_label,
                                      vertex.size=vertex_size,
                                      vertex.color=vertex_color,
                                      vertex.shape=vertex_shape,
                                      edge.arrow.size=0.4, edge.curved=0.1,
                                      layout=layout)
  }
  
  return(graph)
  
  
}