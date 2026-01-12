# FUNCTIONS
library(ggthemr)
library(dendextend)
library(dplyr)
library(pvclust)

#--- Stacked Barplot
custom_colors <- c(
  "#003333", "#006666", "#009999", "#00cccc",
  "#cc3000", "#ff3300", "#ff6600", "#ff9966",
  "#336600", "#669933", "#99cc66", 
  "#333333", "#666666", "#999999", "#cccccc",
  "#663399", "#9966cc", "#cc99ff", "#d8b3ff",
  "#cc9900", "#ff9900", "#ffcc00", "#ffff00",
  "#000099", "#0000ff", "#6699ff", "#66ccff",
  "#660000", "#990000", "#cc0000", "#ff0000",
  "#336666", "#669999", "#99cccc", "#ccffff",
  "#663333", "#996666", "#cc9999", "#ffcccc"
)

custom_colors <- c(
  "#006666", "#009999", "#00cccc",  
  "#cc3000", "#ff3300", "#ff6600", "#ff9966",
  "#336600", "#669933", "#99cc66", 
  "#666666", "#999999", "#cccccc", "#666666", "#999999", "#cccccc", #HVC
  "#663399", "#9966cc", "#cc99ff", "#d8b3ff",
  "#cc9900", "#ff9900", "#ffcc00", "#ffff00",
  "#000099", "#0000ff", "#6699ff", 
  "#660000", "#990000", "#cc0000", "#ff0000",
  "#336666", "#669999", "#99cccc",
  "#663333", "#996666", "#cc9999", "#ffcccc", "#ffccaa"
)

custom_colors <- c(
  colorRampPalette(c("#66cccc", "#99ffff"))(1),
  colorRampPalette(c("#ff9966", "#ffcc99"))(1),
  colorRampPalette(c("#99cc66", "#ccff99"))(1),
  colorRampPalette(c("#999999", "#e6e6e6"))(1),
  colorRampPalette(c("#b399cc", "#e6ccff"))(1),
  colorRampPalette(c("#ffcc66", "#ffff99"))(1),
  colorRampPalette(c("#99ccff", "#cce6ff"))(1),
  colorRampPalette(c("#cc6666", "#ff9999"))(1),
  colorRampPalette(c("#99cccc", "#cce6e6"))(1),
  colorRampPalette(c("#ffcccc"))(1)
)

custom_colors <- c(
  colorRampPalette(c("#006666", "#00cccc"))(3),
  colorRampPalette(c("#cc3000", "#ff9966"))(4),
  colorRampPalette(c("#336600", "#99cc66"))(3),
  colorRampPalette(c("#333333", "#cccccc"))(5),
  colorRampPalette(c("#663399", "#d8b3ff"))(4),
  colorRampPalette(c("#cc9900", "#ffff00"))(4),
  colorRampPalette(c("#000099", "#6699ff"))(3),
  colorRampPalette(c("#660000", "#ff0000"))(4),
  colorRampPalette(c("#336666", "#99cccc"))(3),
  colorRampPalette(c("#663333", "#ffcccc"))(5)
)

custom_colors <- c(
  colorRampPalette(c("#663399", "#d8b3ff"))(4),
  colorRampPalette(c("#336666", "#99cccc"))(3),
  colorRampPalette(c("#336600", "#99cc66"))(3),
  colorRampPalette(c("#cc9900", "#ffff00"))(4),
  colorRampPalette(c("#000099", "#6699ff"))(3),
  colorRampPalette(c("#333333", "#cccccc"))(5),
  colorRampPalette(c("#660000", "#ff0000"))(4),
  colorRampPalette(c("#006666", "#00cccc"))(3),
  colorRampPalette(c("#663333", "#ffcccc"))(5),
  colorRampPalette(c("#cc3000", "#ff9966"))(4)
)


# custom_colors <- c(
#   colorRampPalette(c("#66cccc", "#99ffff"))(1),
#   colorRampPalette(c("#ff9966", "#ffcc99"))(1),
#   colorRampPalette(c("#99cc66", "#ccff99"))(1),
#   colorRampPalette(c("#999999", "#e6e6e6"))(1),
#   colorRampPalette(c("#b399cc", "#e6ccff"))(1),
#   colorRampPalette(c("#ffcc66", "#ffff99"))(1),
#   colorRampPalette(c("#99ccff", "#cce6ff"))(1),
#   colorRampPalette(c("#cc6666", "#ff9999"))(1),
#   colorRampPalette(c("#99cccc", "#cce6e6"))(1),
#   colorRampPalette(c("#ffcccc"))(1)
# )

stackedbarplot_colors = function(meta, groupx, groupfill, fn, horizontal = F, bold = F,
                                 custom_colors = custom_colors){
  
  clkeys = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_keys %>% as.data.frame
  clkeys$size = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_size 
  
  # Add all cell numbers
  spkeys = meta %>% group_by(meta[,groupfill]) %>% group_keys %>% as.data.frame
  colnames(spkeys) = groupfill
  spkeys$size = meta %>% group_by(meta[, groupfill]) %>% group_size
  spkeys = spkeys %>% add_column(tmp = 'AllCells')
  colnames(spkeys)[3] = groupx
  clkeys = rbind(clkeys, spkeys)
  
  colnames(clkeys) = c('cluster', 'variable', 'value')
  
  plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
    geom_bar(position = "fill", stat = "identity") +
    xlab("") +
    ylab("Percentage") +
    theme_classic() +
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme(text=element_text(size=30), axis.text.y = element_text(face = 'bold')) +
    rotate_x_text(45) + coord_flip()
  
  if(horizontal == T){
    plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
      geom_bar(position = "fill", stat = "identity") +
      xlab("") +
      ylab("Percentage") +
      theme_classic() +
      scale_fill_manual(values = custom_colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(text=element_text(size=20)) +
      rotate_x_text(60)
  }
  
  # pdf(paste0(fn, '.pdf'), width = 10, height = 10)
  #print(plt)
  # dev.off()
  
  return(plt)
}

# stackedbarplot_colors = function(meta, groupx, groupfill, fn = NULL, horizontal = F, bold = F,
#                                  custom_colors = custom_colors, level_order = NULL){
#   
#   clkeys = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_keys %>% as.data.frame
#   clkeys$size = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_size 
#   
#   # Add all cell numbers
#   spkeys = meta %>% group_by(meta[,groupfill]) %>% group_keys %>% as.data.frame
#   colnames(spkeys) = groupfill
#   spkeys$size = meta %>% group_by(meta[, groupfill]) %>% group_size
#   spkeys = spkeys %>% add_column(tmp = 'AllCells')
#   colnames(spkeys)[3] = groupx
#   clkeys = rbind(clkeys, spkeys)
#   
#   colnames(clkeys) = c('cluster', 'variable', 'value')
#   
#   # Reapply factor order here
#   if(!is.null(level_order)){
#     clkeys$variable <- factor(clkeys$variable, levels = level_order)
#   }
# 
#   plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) +
#     geom_bar(position = "fill", stat = "identity") +
#     xlab("") +
#     ylab("Percentage") +
#     theme_classic() +
#     scale_fill_manual(values = custom_colors) +
#     scale_y_continuous(labels = scales::percent_format()) +
#     theme(text=element_text(size=30), axis.text.y = element_text(face = 'bold')) +
#     rotate_x_text(45) + coord_flip()
#   
#   if(horizontal == T){
#     plt = plt + coord_flip()
#   }
#   
#   return(plt)
# }

stackedbarplot_colors = function(meta, groupx, groupfill, fn = NULL, horizontal = F, bold = F,
                                 custom_colors = custom_colors, level_order = NULL){
  
  clkeys = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_keys %>% as.data.frame
  clkeys$size = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_size 
  
  # Add all cell numbers
  spkeys = meta %>% group_by(meta[,groupfill]) %>% group_keys %>% as.data.frame
  colnames(spkeys) = groupfill
  spkeys$size = meta %>% group_by(meta[, groupfill]) %>% group_size
  spkeys = spkeys %>% add_column(tmp = 'AllCells')
  colnames(spkeys)[3] = groupx
  clkeys = rbind(clkeys, spkeys)
  
  colnames(clkeys) = c('cluster', 'variable', 'value')
  
  # force factor levels if provided
  if(!is.null(level_order)){
    clkeys$variable <- factor(clkeys$variable, levels = level_order)
  }
  
  plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
    geom_bar(position = "fill", stat = "identity") +
    xlab("") +
    ylab("Percentage") +
    theme_classic() +
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme(text=element_text(size=15), axis.text.y = element_text(face = 'bold'),
          axis.line = element_line(color = "#020202")) +
    rotate_x_text(45) + coord_flip()

  if(horizontal == T){
    plt = plt + coord_flip()
  }
  
  return(plt)
}




stackedbarplot = function(meta, groupx, groupfill, fn, horizontal = F, bold = F){
  
  clkeys = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_keys %>% as.data.frame
  clkeys$size = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_size 
  
  # Add all cell numbers
  spkeys = meta %>% group_by(meta[,groupfill]) %>% group_keys %>% as.data.frame
  colnames(spkeys) = groupfill
  spkeys$size = meta %>% group_by(meta[, groupfill]) %>% group_size
  spkeys = spkeys %>% add_column(tmp = 'AllCells')
  colnames(spkeys)[3] = groupx
  clkeys = rbind(clkeys, spkeys)
  
  colnames(clkeys) = c('cluster', 'variable', 'value')
  
  plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
    geom_bar(position = "fill", stat = "identity") +
    xlab("") +
    ylab("Percentage") +
    theme_classic() +
    # scale_fill_manual(values = custom_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme(text=element_text(size=30), axis.text.y = element_text(face = 'bold')) +
    rotate_x_text(45) + coord_flip()
  
  if(horizontal == T){
    plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
      geom_bar(position = "fill", stat = "identity") +
      xlab("") +
      ylab("Percentage") +
      theme_classic() +
      # scale_fill_manual(values = custom_colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(text=element_text(size=20)) +
      rotate_x_text(60)
  }
  
  # pdf(paste0(fn, '.pdf'), width = 10, height = 10)
  #print(plt)
  # dev.off()
  
  return(plt)
}








