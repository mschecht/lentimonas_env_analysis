#!/usr/bin/env Rscript

# Check if basic packages are installed 
# if not here then install them
#-----------------------------------

packages <- c("tidyverse", "data.table", "parallel", "RSQLite", "rworldmap", "knitr", "phyloseq", "microbiome")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
#---------------
suppressMessages(lapply(packages, library, character.only = TRUE))


# load data
#----------
load("data/processed/mitags.RData")

load("data/interim/mitag_rel_abun.RData")

# Projects contextual data database
db <- RSQLite::dbConnect(drv = SQLite(), dbname = "data/external/contextual_data.db")

dbListTables(db) # list contents of database

# extract TARA Oceans contextual data from database
contex <- dbReadTable(db, "contextual_data_all") %>% as_tibble()
contex_TARA <- dbReadTable(db, "TARA_contex") %>% as_tibble()
contex_OSD <- dbReadTable(db, "OSD_contex") %>% as_tibble()

dbDisconnect(db) # Disconnect from db

# Plot
#----------
prep_map <- function(rel_abun, phylo_lvl_name){
  
  
  ### Test Variables
  #  rel_abun = OSD_mitags_rel_abun$Genus
  #  phylo_lvl_name = "Lentimonas"
  ###
  
  # Convert back to long form for ggplot
  rnames <- rel_abun %>% row.names() # get row names
  
  rel_abun.df <- rel_abun %>% as.data.frame() # convert to dataframe
  
  rownames(rel_abun.df) <- rnames # add in row names (taxa) to dataframe
  
  rel_abun.df %>% 
    as_tibble() %>%
    rownames_to_column(var = "ID") %>% 
    gather(key = "sample_ID", value = "rel_abun", -ID) %>% 
    filter(ID == paste(phylo_lvl_name)) %>%
    left_join(contex %>% select(sample_ID, latitude, longitude, ecoregion)) # Join with metadata
}

Lentimonas_rel_abun_OSD <- prep_map(rel_abun = OSD_mitags_rel_abun$Genus, phylo_lvl_name = "Lentimonas")
Puniceicoccaceae_rel_abun_OSD <- prep_map(rel_abun = OSD_mitags_rel_abun$Family, phylo_lvl_name = "Puniceicoccaceae")
Opitutales_rel_abun_OSD <- prep_map(rel_abun = OSD_mitags_rel_abun$Order, phylo_lvl_name = "Opitutales")
Verrucomicrobiae_rel_abun_OSD <- prep_map(rel_abun = OSD_mitags_rel_abun$Class, phylo_lvl_name = "Verrucomicrobiae")

Lentimonas_rel_abun_TARA <- prep_map(rel_abun = TARA_mitags_rel_abun$Genus, phylo_lvl_name = "Lentimonas")
Puniceicoccaceae_rel_abun_TARA <- prep_map(rel_abun = TARA_mitags_rel_abun$Family, phylo_lvl_name = "Puniceicoccaceae")
Opitutales_rel_abun_TARA <- prep_map(rel_abun = TARA_mitags_rel_abun$Order, phylo_lvl_name = "Opitutales")
Verrucomicrobiae_rel_abun_TARA <- prep_map(rel_abun = TARA_mitags_rel_abun$Class, phylo_lvl_name = "Verrucomicrobiae")



## 4. Visualize relative proportions of `Verrucomicrobiae Opitutales Puniceicoccaceae Lentimonas` at its different phylogenetic levels on the world map


theme_map <- function(...) {
  theme_minimal() +
    theme(
      #text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      ...
    )
}

plot_on_world <- function(df, title){
  
  #df = Lentimonas_rel_abun
  #title = "Relative abundance of Verrucomicrobia Lentimonas in TARA Oceans and OSD"
  # Plot
  map.world <- map_data(map="world")
  
  p <- ggplot() +
    geom_map(data=map.world, map=map.world, aes(map_id=region, x=long, y=lat)) +
    geom_point(data = df, 
               shape =21,
               size = 2.2, 
               alpha = 0.7,
               color = "black",
               aes(x=longitude, y=latitude, fill = rel_abun)) +
    scale_fill_gradient2(low = "#40e0d0", mid = "#ff8c00", high="#ff0080", 
                         #low = "#fcb045", mid = "#fd1d1d", high="#833ab4",
                         trans = "sqrt", 
                         midpoint=1.8,
                         name = "Proportion (%)") +
    ggtitle(paste(title)) +
    xlab("") +
    ylab("") +
    coord_equal(xlim = c(-160, 170)) +
    theme_map(legend.position = "top",
              legend.key.height = unit(0.2, "cm"))
  
  return(p)
}

p_Lentimonas_OSD <- Lentimonas_rel_abun_OSD %>% 
  mutate(rel_abun = (rel_abun) * 100) %>%
  plot_on_world(title = "Relative proportions of Verrucomicrobia Lentimonas in OSD")

p_Puniceicoccaceae_OSD <- Puniceicoccaceae_rel_abun_OSD %>% 
  plot_on_world(title = "Relative proportions of Puniceicoccaceae in OSD")

p_Opitutales_OSD <- Opitutales_rel_abun_OSD %>% 
  plot_on_world(title = "Relative proportions of Opitutales in OSD")

p_Verrucomicrobiae_OSD <- Verrucomicrobiae_rel_abun_OSD %>% 
  plot_on_world(title = "Relative proportions of Verrucomicrobiae Lentimonas in OSD")

plot_on_world1 <- function(df, title) {
  
  #df = Lentimonas_rel_abun
  #title = "Relative abundance of Verrucomicrobia Lentimonas in TARA Oceans and OSD"
  # Plot
  map.world <- map_data(map="world")
  
  p <- ggplot() +
    geom_map(data=map.world, map=map.world, aes(map_id=region, x=long, y=lat)) +
    geom_point(data = df, 
               shape =21,
               size = 2.2, 
               alpha = 0.7,
               color = "black",
               aes(x=longitude, y=latitude, fill = rel_abun)) +
    scale_fill_gradient2(low = "#40e0d0", mid = "#ff8c00", high="#ff0080",
                         #low = "blue", high = "red",
                         trans = "sqrt", 
                         midpoint=1,
                         name = "Proportion (%)") +
    ggtitle(paste(title)) +
    xlab("") +
    ylab("") +
    coord_equal(xlim = c(-160, 170)) +
    theme_map(legend.position = "top",
              legend.key.height = unit(0.2, "cm"))
  
  return(p)
}

p_Lentimonas_TARA <- Lentimonas_rel_abun_TARA %>% mutate(rel_abun = (rel_abun) * 100) %>% filter(rel_abun > 0) %>% 
  plot_on_world1(title = "Relative proportions of Verrucomicrobia Lentimonas in TARA")

p_Puniceicoccaceae_TARA <- Puniceicoccaceae_rel_abun_TARA %>% 
  plot_on_world1(title = "Relative proportions of Puniceicoccaceae in TARA")

p_Opitutales_TARA <- Opitutales_rel_abun_TARA %>% 
  plot_on_world1(title = "Relative proportions of Opitutales in TARA")

p_Verrucomicrobiae_TARA <- Verrucomicrobiae_rel_abun_TARA %>% 
  plot_on_world1(title = "Relative proportions of Verrucomicrobiae Lentimonas in TARA")

ggsave(filename = "reports/figures/p_Lentimonas_OSD.pdf", plot = p_Lentimonas_OSD)

cat("\nPlot 1 done...\n\n")


## 5. Questions about the data

plot_on_world2 <- function(df, title){

#Test variables 
#  df = Lentimonas_rel_abun_OSD
#  title = "Relative abundance of Verrucomicrobia Lentimonas in TARA Oceans and OSD"
###

# Plot
map.world <- map_data(map="world")

p <- ggplot() +
geom_map(data=map.world, map=map.world, aes(map_id=region, x=long, y=lat)) +
geom_point(data = df, 
shape =21,
size = 2.2, 
alpha = 0.7,
color = "black",
aes(x=longitude, y=latitude, fill = "#40e0d0")) +
geom_text(label=df$ecoregion, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T) +
ggtitle(paste(title)) +
xlab("") +
ylab("") +
coord_equal(xlim = c(-160, 170)) +
theme_map(legend.position = "",
legend.key.height = unit(0.2, "cm"))

return(p)
}

# Top 10 Lentimonas
# OSD
Lentimonas_rel_abun_OSD %>% 
ungroup() %>% 
arrange(desc(rel_abun)) %>% 
slice(1:10) %>%
#filter(grepl("OSD113", sample_ID)) %>%
plot_on_world2(title = "Relative proportion of Verrucomicrobia Lentimonas in OSD (top 10 sites)") %>%
ggsave(filename = "reports/figures/osd_top10.pdf")

cat("\nPlot 2 done...\n\n")


# TARA
Lentimonas_rel_abun_TARA %>% 
ungroup() %>% 
arrange(desc(rel_abun)) %>% 
slice(1:0) %>%
plot_on_world2(title = "Relative proportion of Verrucomicrobia Lentimonas in TARA Oceans (top 10 sites)")%>%
ggsave(filename = "reports/figures/tara_top10.pdf")

cat("\nPlot 3 done...\n\n")


# OSD
# How much?
Lentimonas_rel_abun_OSD_max_prop <- (Lentimonas_rel_abun_OSD %>% 
ungroup() %>% 
arrange(desc(rel_abun)) %>% 
slice(1:10) %>%
.$rel_abun) * 100

osd_high_sam <- Lentimonas_rel_abun_OSD%>% 
arrange(desc(rel_abun)) %>% 
slice(1:10) %>% .$sample_ID


## Heat map of various Genus of Lentimonas in samples

# Sample with top 10 most Puniceicoccaceae
OSD_samples_lentimonas_10 <- Puniceicoccaceae_rel_abun_OSD %>% 
ungroup() %>% 
arrange(desc(rel_abun)) %>% 
slice(1:10) %>%
.$sample_ID

TARA_samples_lentimonas_10 <- Puniceicoccaceae_rel_abun_TARA %>% 
ungroup() %>% 
arrange(desc(rel_abun)) %>% 
slice(1:10) %>%
.$sample_ID

plot_cool_heatmap<- function(samples, mitags){

#Test variables
###
#  samples = OSD_samples_lentimonas_10
#  mitags = osd_mitags
###

# ggplot stuff
base_breaks <- function(n = 10){
function(x) {
axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
}
}

textcol <- "grey40"

# Plot samples with top 10 proportions of Lentimonas
mitags %>% 
filter(Family == "Puniceicoccaceae") %>%
filter(label %in% samples) %>%
group_by(label, Genus) %>%      
add_count() %>%                 
ungroup() %>% 
group_by(label) %>%                             
mutate("miTAG_count_per_sample" = sum(n)) %>%   
mutate(rel_abun = n/miTAG_count_per_sample) %>%
ungroup() %>%

ggplot(aes(y=label,x=Genus, fill = rel_abun))+ ##### weird that
geom_tile()+
#redrawing tiles to remove cross lines from legend
geom_tile(colour="white",size=0.25, show.legend = FALSE)+
#remove axis labels, add title
labs(x="",y="",title="")+
#remove extra space
scale_y_discrete(expand=c(0,0))+
#custom breaks on x-axis
scale_x_discrete(expand=c(0,0))+
#custom colours for cut levels and na values
#scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
# "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
viridis::scale_fill_viridis(option = "D", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#mark year of vaccination
#geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
#equal aspect ratio x and y axis
coord_fixed() +
#set base size for all font elements
theme_grey(base_size=10)+
#theme options
theme(
#legend.position = "bottom",
#axis.text.x = element_text(angle = 90, hjust = 1),
#remove legend title
#legend.title=element_blank(),
legend.title=element_text("test"),
#remove legend margin
legend.spacing = grid::unit(0,"cm"),
#change legend text properties
legend.text=element_text(colour=textcol,size=7,face="bold"),
#change legend key height
legend.key.height=grid::unit(0.8,"cm"),
#set a slim legend
legend.key.width=grid::unit(0.8,"cm"),
#set x axis text size and colour
axis.text.x=element_text(hjust = 1, vjust = 0.5, colour=textcol, angle = 90, size = 6),
axis.ticks.y = element_blank(),
#axis.text.y = element_blank(),
#set y axis text colour and adjust vertical justification
axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
#change axis ticks thickness
axis.ticks=element_line(size=0.4),
#change title font, size, colour and justification
plot.title=element_blank(),
#remove plot background
plot.background=element_blank(),
#remove plot border
panel.border=element_blank())
}

p_tara_heat<- plot_cool_heatmap(TARA_samples_lentimonas_10, tara_mitags)

p_osd_heat <- plot_cool_heatmap(OSD_samples_lentimonas_10, osd_mitags)

ggsave(filename = "reports/figures/p_osd_letimonas_heat.pdf", plot = p_osd_heat)

cat("\nPlot 4 done...\n\n")
