require(tidyverse)
require(pbapply)
require(parallel)
require(lisi)
require(emdist)

# require(pbmcapply)

set.seed(222)

# Function to compute lisi scores from shuffled data
shuffled_lisi <- function(meta = Harmony_Meta, embeds = HarmonyUMAP_embeds, shuf_num = 1000){
  # Shuffle Meta brain regions within individual cell types
  shuffled_lists <- replicate(shuf_num, .shuffle_within_types(meta), simplify = FALSE)
  
  # Combine shuffled brain regions lists into a single data frame
  df_BR_shuffled <- bind_cols(sapply(shuffled_lists, `[`, "BR_shuf"), .id = NULL)
  
  # Combine Harmony Meta with shuffled data frame
  Harmony_Meta_shuffled <- cbind(meta, df_BR_shuffled)
  
  # Extract column names to compute LISI scores for
  column_names <- colnames(Harmony_Meta_shuffled)[-(1:2)]
  
  # compute lisi socres in parallel
  Region_shuffled_lisi <- pblapply(column_names, \(col) {
    compute_lisi(embeds, Harmony_Meta_shuffled, col)
  }, cl = detectCores())|>
    as.data.frame()|>
    rename_with(~ paste0(., "_lisi"))
  
  return(list(Meta_shuffled = Harmony_Meta_shuffled,
              lisi_shuffled = Region_shuffled_lisi
  )
  )
}

# Function to reshape lisiScores data frame for plotting 
reshape_lisi <- function(df = lisiScores, cell_type = "all") {
  df <- df |>
    select(-c(Region, RegionAcronym))
  
  # Determine cell types to process:
  v <- unique(df$CellSubclass)
  if (identical(cell_type, "all")) {
    val <- v
  } else if (all(cell_type %in% v)) {
    val <- cell_type
  } else {
    cat ("cell_type has to be either 'all' or a list of valid cell types:", "\n", v, "\n")
    return(NULL)
  }
  
  # Reshape based on selected cell types:
  lisiScores_reshaped <- pblapply(val, \(val) {
    .process_group(df, val)
  }, cl = detectCores()) |>
    bind_rows()
  
  return(lisiScores_reshaped)
}

# Function to plot density plot
plot_lisi <- function(df, input = c("all", "original")[1], save_build = FALSE,
                      print_plot = TRUE, save_plot = FALSE, 
                      saveto = "figs/lisi_scores.pdf"){
  if (identical(input, "original")){
    df <- df |>
      filter(source == "original" )
    plot <- ggplot(df, aes(lisi, fill=CellSubclass, color=CellSubclass)) + 
      geom_density() + 
      facet_grid(CellSubclass~., scales="free") + 
      # scale_color_manual(name = "", values = text_colors) +
      # scale_fill_manual(name = "", values = text_colors) +
      theme(legend.position="none",
            strip.background = element_blank(),
            strip.text.y = element_text(angle=0, hjust=0)) + 
      # xlim(.8,2.2) + 
      # scale_x_continuous(breaks=c(1,2), limits=c(.8,2.2)) + 
      scale_y_continuous(limits=c(0,0.5)) + 
      labs(x="Region LISI, UMAP Space")
  }else{
    plot <- ggplot(df, aes(lisi, fill = source, color = source)) +
      geom_density(alpha = 0.5, position = "identity") +
      facet_grid(CellSubclass ~ ., scales = "free") +  # Use facet_grid for vertical facets
      theme(
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0),
        legend.position = "top",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0)
      ) +
      scale_color_manual(values = c("original" = "blue", "shuffle" = "gray")) +
      scale_fill_manual(values = c("original" = "blue", "shuffle" = "gray")) +
      labs(x = "Region LISI, UMAP Space", y = "Density") +
      guides(fill = guide_legend(title = "Source"), color = guide_legend(title = "Source"))
  }
  
  if (save_plot) {
    ggsave(saveto, plot, device = "pdf")
  }
  
  if (print_plot) {
    print(plot)
  }
  
  if (save_build) {
    return (ggplot_build(plot)$data[[1]])
  }
}

# Function to retrieve kernel density estimate from ggplot
KDE_lisi <- function(df = p){
  df <- p|>
    select(PANEL, group, x , y)|>
    rename(CellSubclass = PANEL, Source = group)
  
  cell_types <- setNames(c("Astrocytes", "Endothelial" , "GABAergic", "Glutamatergic",
                           "Microglia", "Oligodendrocytes", "OPCs"), 1:7)
  
  df$CellSubclass <- cell_types[df$CellSubclass]
  
  df$Source <- replace(df$Source, df$Source == 1, "original")
  df$Source <- replace(df$Source, df$Source == 2, "shuffle")
  
  return(df)
}

# Function to calculate earth mover's distance(EMD)/ dissimilarity across different CellSubclass 
emd_lisi <- function(df = KDE){
  cell_types <- unique(df$CellSubclass)
  emd_res <- pbsapply(cell_types, \(val) {
    original_xy <- df |> filter(CellSubclass == val, Source == "original") |> select(3:4) |> as.matrix()
    shuffle_xy <- df |> filter(CellSubclass == val, Source == "shuffle") |> select(3:4) |> as.matrix()
    emd <- emd(original_xy, shuffle_xy)
  }) |> 
    cbind() |>
    as.data.frame() |>
    setNames("EMD") |>
    arrange(desc(EMD))
}

# Function to shuffle brain regions within each cell type
.shuffle_within_types <- function(df) {
  df |>
    group_by(CellSubclass) |>
    mutate(BR_shuf = sample(RegionAcronym)) |>
    ungroup()
}

# Function to reshape data frame for each cell type
.process_group <- function(df, val){
  df_reshaped <-  df|>
    filter(CellSubclass == val)|>
    pivot_longer(cols = -CellSubclass, names_to = "source", values_to = "lisi") |>
    mutate(source = if_else(str_starts(source, "BR_shuf"), "shuffle", "original"))|>
    arrange(CellSubclass, source)
}


#################################################################

# to be optimized 

# Function to calculate statistics significance for cell types
stats_lisi <- function(df = lisiScores, cell_type = "all") {
  df <- df |>
    select(-c(Region, RegionAcronym)) 
  
  # Determine cell types to process:
  v <- unique(df$CellSubclass)
  if (identical(cell_type, "all")) {
    val <- v
  } else if (all(cell_type %in% v)) {
    val <- cell_type
  } else {
    cat ("cell_type has to be either 'all' or a list of valid cell types:", "\n", v, "\n")
    return(NULL)
  }
  
  # calculate column Means
  nCores = detectCores()
  stats <- pblapply(val, \(val) {
    .colMenas_group(df, val)
  }, cl = nCores)|>
    bind_rows()
  rownames(stats) <- val
  
  # calculate P values
  i <- nrow(stats)
  p_value <- sapply(1:i, \(i){ 
    mean(stats[i,-1] < as.vector(stats[i,1]))
  })
  
  return(data.frame (CellSubclass = val,
                     p_value = p_value)
  )
}

# Function to calculate column mean for each cell type
.colMenas_group <- function(df, val){
  stat <-  df |>
    filter(CellSubclass == val)|>
    select(-CellSubclass)|>
    pbapply(2, mean, cl = nCores)
}    



###--- earthmovers distance
KDE_lisi <- function(df = p){
  df <- df|>
    select(PANEL, group, x , y)|>
    rename(CellSubclass = PANEL, Source = group)
  
  cellClassNames <- levels(factor(combined_lisiScores$CellSubclass))
  cell_types <- setNames(cellClassNames, 1:length(cellClassNames))
  
  df$CellSubclass <- cell_types[df$CellSubclass]
  
  df$Source <- replace(df$Source, df$Source == 1, "nonSong")
  df$Source <- replace(df$Source, df$Source == 2, "song")
  
  return(df)
}


emd_lisi <- function(df = KDE){
  cell_types <- unique(df$CellSubclass)
  emd_res <- pbsapply(cell_types, \(val) {
    original_xy <- df |> filter(CellSubclass == val, Source == "nonSong") |> select(3:4) |> as.matrix()
    shuffle_xy <- df |> filter(CellSubclass == val, Source == "song") |> select(3:4) |> as.matrix()
    emd <- emd(original_xy, shuffle_xy)
  }) |> 
    cbind() |>
    as.data.frame() |>
    setNames("EMD") |>
    arrange(desc(EMD))
}


