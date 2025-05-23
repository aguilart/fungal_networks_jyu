
rm(list = ls()); gc()

# loading packages
library(tidyverse)
library(readxl)
library(tidytext)
library(vegan)
library(ranger)

# loading functions
my_theme <-  theme_minimal(base_size = 18) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.position = "right",
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dotted', color = 'gray'),
    panel.grid.minor = element_blank()
  )

imp <- function(ord, n){
  round(summary(ord)[["cont"]][["importance"]][2,n],2)*100
}

merge_rda <- function(rda_, dat){
  m <-
    sum(
      length(rda_$CCA$eig),
      length(rda_$CA$eig))
  
  test3 <- as.data.frame(scores(rda_,
                                display = "sites",
                                scaling = "sites",
                                choices=c(1:m)))
  
  test3 <- rownames_to_column(test3) %>% rename(name_col = rowname)
  
  
  merged <- left_join(dat,
                      test3, by = "name_col")
  merged}


get_arrows <-
  function(datos, n_traits, trait_codes = F) {
    m <- sum(
      length(datos$CCA$eig),
      length(datos$CA$eig))
    contributions <- as.data.frame(
      scores(datos, display = "species",
             choices = c(1:m),
             scaling =  0));
    contributions$traits <- rownames(contributions)
    traits <- unique(c(contributions$traits[order(contributions[,1]^2, decreasing = T)][1:n_traits],
                       contributions$traits[order(contributions[,2]^2, decreasing = T)][1:n_traits]))
    contributions <- contributions[which(contributions$traits%in%traits),]
    if(trait_codes == T){
      contributions <- left_join(contributions, id_mapping, by = "traits")
    }else{contributions}
  }

# work in progress
# get_arrows <-
#   function(datos, n_traits, trait_codes = F, x, y) {
#     m <- sum(
#       length(datos$CCA$eig),
#       length(datos$CA$eig))
#     contributions <- as.data.frame(
#       scores(datos, display = "species",
#              choices = c(1:m),
#              scaling =  0));
#     contributions$traits <- rownames(contributions)
#     traits <- unique(c(contributions$traits[order(contributions[,x]^2, decreasing = T)][1:n_traits],
#                        contributions$traits[order(contributions[,y]^2, decreasing = T)][1:n_traits]))
#     contributions <- contributions[which(contributions$traits%in%traits),]
#     if(trait_codes == T){
#       contributions <- left_join(contributions, id_mapping, by = "traits")
#     }else{contributions}
#   }




# specifying folder
data_location2 <- "C:/Users/mafekete/Dropbox/fungi/results_after_clearing/"

#loading data
all_data_scaled_a <- readRDS(paste(data_location2, "all_data_scaled_a.RDS", sep = ""))

# create a template for meta data entry
write.csv(data.frame(name_col = all_data_scaled_a$name_col),
          paste(data_location2,"meta_data.csv", sep = ""), row.names = F)

# After filling the csv above, call it

meta_data <- read.csv(paste(data_location2, "meta_data_added.csv", sep = ""), sep = ";")

# Merge data with metadata

all_data_scaled_a <- left_join(all_data_scaled_a, meta_data)

# Specifiy the variables to use for the PCA
scaled_variables <- c(#"Hyphal_main_width_median",
  
  "Hyphal_main_width_mean_of_MST", # To include only the main roads # 1
  "Hyphal_main_length_median_of_MST", # 2
  "Hyphal_sec_length_median", # 3
  "Hyphal_sec_width_mean",   # 4
  "Hyphal_angle", # 5
  
  "Hyphal_density_main", # 6
  
  "colony_Circularity", # 7
  
  
  # Measures of heterogeneity of the main roads
  "accessibility_skewness", # 8
  "length_skewness_main_of_MST", # 9
  "width_skewness_main_of_MST", # 10
  
  
  # Transport efficiency:
  "Reff_tip",  # 11
  #"Geff_MST",             
  "Geff_MST_no_c", # 12   
  
  # Loopiness index
  "alpha_coeff", # 13
  
  "Mycelia_Vol_MST", # 14
  
  # Robustness measures:
  
  "spatial_random_AUC_mean", # 15
  "spatial_random_AUC_cv" , # 16
  
  "AUC_Length_asc", # 17                   
  "AUC_Length_des", # 18
  "AUC_Resistance_asc", # 19               
  
  "AUC_Width_asc",  # 20
  "AUC_Width_des"#, # 21
)

# 


library(RColorBrewer)
custom_colors <- brewer.pal(n = 7, name = "Set2")



#exp1 <- grep("ctrl|coll|wl1|nema|wl1|wl2|mill", all_data_scaled_a$Grazer)

datos1 <- as_tibble(all_data_scaled_a[, c("name_col", scaled_variables)]) %>%
  column_to_rownames(var = "name_col") #%>% 
  #filter(!grepl("ceuthospora_3|ceuthospora_4",all_data_scaled_a$name_col))

pca_exp1 <-
  rda(datos1[, scaled_variables],
      scale = TRUE, data = datos1)

pca_plot1 <-
  merge_rda(pca_exp1,
            all_data_scaled_a[, c("name_col", "treatment", scaled_variables)] ) %>% 
              #filter(!grepl("ceuthospora_3|ceuthospora_4",all_data_scaled_a$name_col))) %>% 
  ggplot()+
  aes(x = PC1, y = PC2) +
  geom_point(aes(colour = treatment), size = 5, alpha = 0.7) +
  geom_text(aes(label = name_col)) +
  #facet_wrap(.~day, nrow = 1) +
  #stat_ellipse(aes(group = species), alpha = 0.2, color = "black", linetype = "solid", size = 1) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  # labs(
  #   x = "PC1",
  #   y = "PC2"
  # ) + 
  labs(x = paste("PC1", imp(pca_exp1,1), "%", sep = " "),
       y = paste("PC2", imp(pca_exp1,2),"%", sep = " ")) +
  my_theme

pca_plot1




# Get arrows to plot on top of the pca

arrows <-
get_arrows(pca_exp1, n_traits = 5)

s <- 0.3

pca_plot1 + 
geom_segment(aes(x=0, y=0, xend= arrows[,1]/s, yend=arrows[,2]/s),
             linewidth = 1, arrow = arrow(),
             data = arrows)+
  geom_text_repel(size=4,aes(x = arrows[,1]/s, y = arrows[,2]/s,
                       label = traits,
                       fontface="bold"),
            hjust = 1.1, vjust =-0.5, data = arrows)


ggsave(paste(data_location2,"Figures//PCA_exp1_day.png", sep = ""), pca_plot1, width = 16, height = 8, dpi = 300)

#############################################################################################################

# Making one big plot with the results from the permanova and betadisper
pca_function <- 
  function(dat, variables){
    
    rda <- rda(dat[,variables] ~ Condition(PC1 + PC2), scale = T, data = dat)
    m <-
      sum(
        length(rda$CCA$eig),
        length(rda$CA$eig))
    
    test3 <- as.data.frame(scores(rda,
                                  display = "sites",
                                  scaling = "sites",
                                  choices=c(1:m)))
    
    test3 <- rownames_to_column(test3) %>% rename(name_col = rowname)
    dat <- rownames_to_column(dat) %>% rename(name_col = rowname)
    
    merged <- left_join(dat[,c("name_col" ,"Species", "Grazer")],
                        test3, by = "name_col")
    merged
    
  }


experiment_end <- all_data_scaled_a[, c("Species", "Grazer", "day", "name_col",
                                        "ID", "PC1", "PC2", scaled_variables)]
experiment_end <- experiment_end %>% filter(day == 8)

step1 <- # calculate the big PCA taking into account initial differences
  pca_function(dat = experiment_end %>% 
                 column_to_rownames(var = "name_col"),
               variables = scaled_variables)

step2 <- step1 %>% # calculating the centroids of each species x grazer combinations
  group_by(Species, Grazer) %>% 
  summarise_at(vars(PC1:PC2), mean) %>% ungroup()

# estimate the dispersion between each grazer centroid (among fungal species)
step3 <- dispersion_function(step2, perm = F)

# Extract the median distances between each fungi in the respective grazer and get as a dataframe for boxplots
step4 <- get_tables(step3)

#### Making boxplots of dissimilarity ####
step4 %>% 
  ggplot() +
  aes(treatments, distances) +
  geom_boxplot(aes(fill = treatments), alpha = 0.5, outliers = F) +
  # geom_jitter(aes(shape = ), size = 3) +
  labs(y = "Morphological disimilarity\namong the four species") +
  scale_fill_manual(
    values = c("a_ctrl" = "darkgray",
               "coll" = "red",
               "mill" = "blue",
               "wl1" = "purple",
               "wl2" = "#D55E00",
               "nema" = "green"
    )) +
  scale_x_discrete(
    labels = c("a_ctrl" = "Control",
               "coll" = "Collembola",
               "mill" = "Millipedes",
               "nema" = "Nematodes",
               "wl1" = "Woodlice (1)",
               "wl2" = "Woodlice (2)") # Custom labels
  ) +
  theme_minimal() +
  theme_minimal(base_size = 18) +
  theme(
    #strip.text = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 20),
    legend.position = "none",
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dotted', color = 'gray'),
    panel.grid.minor = element_blank()
  )
#####

# The boxplots above make a really good point, however I think that readers would rather see the actual dispersion
# along the axes. So an easy solution that is coherent with the analysis above is to simply plot the big PCA and 
# show the amount of variation

centroids <- step2 %>%
  group_by(Grazer) %>%
  summarise(across(PC1:PC2, mean, .names = "centroid_{.col}"))

facet_labels = c("a_ctrl" = "Control",
                 "coll" = "Collembola",
                 "mill" = "Millipedes",
                 "nema" = "Nematodes",
                 "wl1" = "Woodlice (1)",
                 "wl2" = "Woodlice (2)") # Custom labels



cowplot::plot_grid(
  step2 %>% 
    ggplot()+
    aes(x = PC1, y = PC2) +
    facet_wrap(.~Grazer, nrow = 3, ncol = 2,labeller = labeller(Grazer = facet_labels)) +
    stat_ellipse(aes(group = Grazer, fill = Grazer), alpha = 0.2,
                 linetype = "solid", size = 1, geom = "polygon",
                 level = 0.7, show.legend = c(fill = F)) +
    geom_point(aes(shape = Species), size = 5) +
    geom_text(aes(centroid_PC1, centroid_PC2, color= Grazer,
                  label = "\u2605"), size = 8, show.legend = c(color = F), data = centroids) +
    # scale_color_manual(values = custom_colors) +
    # scale_fill_manual(values = custom_colors) +
    labs(
      x = "PC1",
      y = "PC2"
    ) + 
    scale_fill_manual(
      values = c("a_ctrl" = "darkgray",
                 "coll" = "red",
                 "mill" = "blue",
                 "wl1" = "purple",
                 "wl2" = "#D55E00",
                 "nema" = "green"
      )) +
    
    scale_color_manual(
      values = c("a_ctrl" = "darkgray",
                 "coll" = "red",
                 "mill" = "blue",
                 "wl1" = "purple",
                 "wl2" = "#D55E00",
                 "nema" = "green"
      )) +
    
    theme_minimal(base_size = 18) +
    theme(
      strip.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 20),
      legend.position = "bottom",
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'dotted', color = 'gray'),
      panel.grid.minor = element_blank()
    ),
  
  table_grazer_effect_size %>% 
    filter(grepl("Grazer", source)) %>%
    filter(p_value < 0.06) %>%
    arrange(desc(R2)) %>% 
    ggplot() +
    aes(x = reorder(id, R2), y = R2) +
    geom_bar(aes(fill = id) ,stat = "identity", alpha = 0.7) +
    #labs(title = "Deviation of treatment to control") +
    ylab(label = "Proportion of explained variance by grazer") +
    
    scale_fill_manual(
      values = c("millipedes" = "blue",
                 "woodlice_1" = "purple",
                 "woodlice_2" = "#D55E00",
                 "nematodes" = "green")) +
    theme_minimal(base_size = 18) +
    
    theme(axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          legend.position = "bottom",
          axis.title.y = element_blank()) +
    coord_flip(),
  ncol = 1, rel_heights = c(1, 0.3)  )

ggsave("Figures//similarities_grazer.pdf", plot = last_plot(),
       width = 12, height = 15, dpi = 300, units = "in")




resultados_1 %>% 
  #filter(grepl("Grazer_x", source)) %>%
  filter(Source != "Total") %>%
  arrange(desc(R2)) %>% 
  ggplot() +
  aes(x = reorder(Source, R2), y = R2) +
  geom_bar(stat = "identity") +
  #labs(title = "Deviation of treatment to control") +
  ylab(label = "Proportion of explained variance of treatment") +
  xlab(label = "Treatment") +
  theme_minimal(base_size = 18) +
  theme(
    #strip.text = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 40, face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 30),
    legend.position = "right",
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dotted', color = 'gray'),
    panel.grid.minor = element_blank()
  ) +  coord_flip()

ggsave("explained_variation.png", last_plot(), width = 16, height = 8, dpi = 300)


all_data_scaled_a$Grazer[which(all_data_scaled_a$Species == "Phallus impudicus"&
                                 grepl("ctrl", all_data_scaled_a$Grazer))] <- "a_ctrl2"

exp2 <- grep("ctrl2|nema|wl2|mill", all_data_scaled_a$Grazer)




# Trajectories using mean PC axes

# With PCA

sps_experiment <- all_data_scaled_a[, c("Species", "Grazer", "day", "name_col", "ID",
                                        "PC1", "PC2", scaled_variables)]

sps_experiment <- sps_experiment %>% filter(day < 9)

pca_function <- 
  function(dat, variables){rda(dat[,variables] ~ Condition(PC1 + PC2), scale = T, data = dat)}

# pca_function <- 
#   function(dat, variables){rda(dat[,variables] , scale = T, data = dat)}

PCS <- c("PC1", "PC2")

day_interactions <-
  lapply(list(
    #Hypholoma
    sps_experiment %>%  filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),
    
    #Phanerochaete
    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),
    
    #Phallus
    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),
    
    #Resinicium
    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col")),
    
    pca_function, variables = scaled_variables)
names(day_interactions ) <-  c("Hypholoma", "Phanerochaete", "Phallus","Resinicium")





# data <- 
# merge_rda(day_interactions$Hypholoma,
#           sps_experiment %>% filter(grepl("mill|wl2|nema|a_ctrl2", Grazer)) %>% filter(grepl("Hypho", Species)) )

#title <-  ggdraw() + draw_label("With RDA Grazer:day")

sps_experiment2 <- sps_experiment
sps_experiment2$Grazer <- gsub("a_ctrl", "z_ctrl", sps_experiment2$Grazer)


resi <- merge_rda(day_interactions$Resinicium,
                  sps_experiment2 %>% select(!all_of(PCS)) %>%
                    filter(grepl("coll|mill|wl2|wl1|nema|z_ctrl", Grazer)) %>%
                    filter(grepl("Resi", Species)) )

hypho <- merge_rda(day_interactions$Hypholoma,
          sps_experiment2 %>% select(!all_of(PCS)) %>%
            filter(grepl("coll|mill|wl2|wl1|nema|z_ctrl", Grazer)) %>%
            filter(grepl("Hypho", Species)) )

phane <- merge_rda(day_interactions$Phanerochaete,
                   sps_experiment2 %>% select(!all_of(PCS)) %>%
                     filter(grepl("coll|mill|wl2|wl1|nema|z_ctrl", Grazer)) %>%
                     filter(grepl("Phane", Species)) )

phall <- merge_rda(day_interactions$Phallus,
                   sps_experiment2 %>% select(!all_of(PCS)) %>%
                     filter(grepl("coll|mill|wl2|wl1|nema|z_ctrl", Grazer)) %>%
                     filter(grepl("Phall", Species)) )

####

effect_resi <- get_effect_size(resi)
effect_hypho <- get_effect_size(hypho)
effect_phane <- get_effect_size(phane)
effect_phall <- get_effect_size(phall)

####

for_resi <- prepare_data(resi)
for_hypho <- prepare_data(hypho)
for_phane <- prepare_data(phane)
for_phall <- prepare_data(phall)

#### 

library(cowplot)
plot_grid(
  
  plot_grid(
  ggdraw() + draw_label("Resinicium bicolor", fontface = 'italic', size = 12, hjust = 0.5, vjust = 1),
  plot_grid(
  plot_pc_arrow_ribbons(
   data = for_resi,
   var = "PC1",
   colores = c(
     "z_ctrl" = "black",  
     "mill" = "lightgray", 
     "nema" = "lightgray",
     "coll" = "lightgray",
     "wl1" = "purple",
     "wl2" = "#D55E00"
   )),  plot_pc_arrow_ribbons(
    data = for_resi,
    var = "PC2",
    colores = c(
      "z_ctrl" = "black",  
      "mill" = "lightgray", 
      "coll" = "lightgray",
      "nema" = "lightgray",
      "wl1" = "purple",
      "wl2" = "lightgray"
    )), ncol = 2), ncol =1, rel_heights = c(0.05, 1)) ,
  
  plot_grid(
    ggdraw() + draw_label("Hypholoma fasciculare", fontface = 'italic', size = 12, hjust = 0.5),
    plot_grid(
  plot_pc_arrow_ribbons(
    data = for_hypho,
    var = "PC1",
    colores = c(
      "z_ctrl" = "black",  
      "mill" = "lightgray", 
      "nema" = "lightgray",
      "coll" = "lightgray",
      "wl1" = "lightgray",
      "wl2" = "lightgray"
    )),  plot_pc_arrow_ribbons(
      data = for_hypho,
      var = "PC2",
      colores = c(
        "z_ctrl" = "black",  
        "mill" = "blue", 
        "nema" = "green",
        "coll" = "lightgray",
        "wl1" = "lightgray",
        "wl2" = "#D55E00"
      )), ncol = 2), ncol =1, rel_heights = c(0.05, 1)) ,
  
  
  plot_grid(
    ggdraw() + draw_label("Phanerochaete velutina", fontface = 'italic', size = 12, hjust = 0.5),
    plot_grid(
  plot_pc_arrow_ribbons(
    data = for_phane,
    var = "PC1",
    colores = c(
      "z_ctrl" = "black",  
      "mill" = "blue", 
      "nema" = "green",
      "coll" = "lightgray",
      "wl1" = "purple",
      "wl2" = "#D55E00"
    )),  plot_pc_arrow_ribbons(
      data = for_phane,
      var = "PC2",
      colores = c(
        "z_ctrl" = "black",  
        "mill" = "blue", 
        "nema" = "lightgray",
        "coll" = "lightgray",
        "wl1" = "lightgray",
        "wl2" = "#D55E00"
      )), ncol = 2), ncol =1, rel_heights = c(0.05, 1)) ,
  
  plot_grid(
    ggdraw() + draw_label("Phallus impudicus", fontface = 'italic', size = 12, hjust = 0.5),
    plot_grid(
  plot_pc_arrow_ribbons(
    data = for_phall,
    var = "PC1",
    colores = c(
      "z_ctrl" = "black",  
      "mill" = "blue", 
      "nema" = "green",
      "coll" = "lightgray",
      "wl1" = "lightgray",
      "wl2" = "#D55E00"
    )),  plot_pc_arrow_ribbons(
      data = for_phall,
      var = "PC2",
      colores = c(
        "z_ctrl" = "black",  
        "mill" = "lightgray", 
        "nema" = "lightgray",
        "coll" = "lightgray",
        "wl1" = "lightgray",
        "wl2" = "lightgray"
      )),ncol = 2), ncol =1, rel_heights = c(0.05, 1)) ,
  
  ggplot() + theme_void() + # Blank plot
    annotate("segment", x = 0.15, xend = 0.25, y = 0.50, yend = 0.50, color = "black", linewidth = 2) +
    annotate("text", x = 0.27, y = 0.50, label = "Control", hjust = 0, size = 3) +
    
    annotate("segment", x = 0.45, xend = 0.55, y = 0.50, yend = 0.50, color = "blue", linewidth = 2) +
    annotate("text", x = 0.57, y = 0.50, label = "Millipedes", hjust = 0, size = 3) +
    
    annotate("segment", x = 0.8, xend = 0.9, y = 0.50, yend = 0.50, color = "purple", linewidth = 2) +
    annotate("text", x = 0.92, y = 0.50, label = "Woodlice (1)", hjust = 0, size = 3) +
    # Second row
    annotate("segment", x = 0.15, xend = 0.25, y = 0.49, yend = 0.49, color = "#D55E00", linewidth = 2) +
    annotate("text", x = 0.27, y = 0.49, label = "Woodlice (2)", hjust = 0, size = 3) +
    
    annotate("segment", x = 0.45, xend = 0.55, y = 0.49, yend = 0.49, color = "green", linewidth = 2) +
    annotate("text", x = 0.57, y = 0.49, label = "Nematodes", hjust = 0, size = 3) +
    
    annotate("segment", x = 0.8, xend = 0.9, y = 0.49, yend = 0.49, color = "red", linewidth = 2) +
    annotate("text", x = 0.92, y = 0.49, label = "Collembola", hjust = 0, size = 3),
  
  ncol = 1, rel_heights = c(1, 1, 1, 1, 0.17))


ggsave("PC_axes_selected_fungi.png", plot = last_plot(),
       width = 7, height = 9, dpi = 300, units = "in")

ggsave("PC_axes_selected_fungi.pdf", plot = last_plot(),
       width = 7, height = 9, dpi = 300, units = "in")





for_resi_nema <- prepare_data(merge_rda(day_interactions$Resinicium,
                                   sps_experiment2 %>% select(!all_of(PCS)) %>%
                                     filter(grepl("nema|z_ctrl", Grazer)) %>%
                                     filter(grepl("Resi", Species)) )
                              
                              )

plot_grid(
  ggdraw() + draw_label("Resinicium bicolor", fontface = 'italic', size = 12, hjust = 0.5, vjust = 1),
  plot_grid(
    plot_pc_arrow_ribbons(
      data = for_resi_nema,
      var = "PC1",
      colores = c(
        "z_ctrl" = "black",  
        "mill" = "lightgray", 
        "nema" = "green",
        "coll" = "lightgray",
        "wl1" = "purple",
        "wl2" = "#D55E00"
      )),  plot_pc_arrow_ribbons(
        data = for_resi_nema,
        var = "PC2",
        colores = c(
          "z_ctrl" = "black",  
          "mill" = "lightgray", 
          "coll" = "lightgray",
          "nema" = "green",
          "wl1" = "purple",
          "wl2" = "#D55E00"
        )), ncol = 2), ncol =1, rel_heights = c(0.05, 1))



# Heatmaps


get_contributions_variables <- function(rda){
  contributions <- as.data.frame(
    scores(rda, display = "species",
           choices = c(1:3),
           scaling = 0));
  
  # names(contributions)[1]<-"Axis1"
  # names(contributions)[2]<-"Axis2"
  
  contributions$PC1 <- 100*(contributions$PC1^2) #change to PC for comparison
  contributions$PC2 <- 100*(contributions$PC2^2) #change to PC for comparison
  contributions$traits <- rownames(contributions)
  rownames(contributions) <- NULL
  contributions
}


get_contributions_exp_variables <- function(rda, dat){
  m <-
    sum(
      length(rda$CCA$eig),
      length(rda$CA$eig))
  
  test3 <- as.data.frame(scores(rda,
                                display = "bp",
                                scaling = 0,
                                choices=c(1:3)))
  
  test3$RDA1 <- 100*(test3$RDA1^2)
  test3$RDA2 <- 100*(test3$RDA2^2)
  test3$RDA3 <- 100*(test3$RDA3^2)
  test3$traits <- rownames(test3)
  rownames(test3) <- NULL
  #test3$traits[3] <- "grazer_x_day"
  test3 
  
}


for_heatmap <- bind_rows(
  get_contributions_variables(species_single_comps$Hypholoma_mill) %>% mutate(comp = "Hypholoma_mill"),
  get_contributions_variables(species_single_comps$Hypholoma_nema) %>% mutate(comp = "Hypholoma_nema"),
  get_contributions_variables(species_single_comps$Hypholoma_wl2) %>% mutate(comp = "Hypholoma_wl2"),
  
  get_contributions_variables(species_single_comps$Phanerochaete_wl1) %>% mutate(comp = "Phanerochaete_wl1"),
  # get_contributions_variables(species_single_comps$Phanerochaete_wl2) %>% mutate(comp = "Phanerochaete_wl2"),
  # get_contributions_variables(species_single_comps$Phanerochaete_nema) %>% mutate(comp = "Phanerochaete_nema"),
  
  get_contributions_variables(species_single_comps$Resinicium_nema) %>% mutate(comp = "Resinicium_nema"),
  get_contributions_variables(species_single_comps$Resinicium_wl2) %>% mutate(comp = "Resinicium_wl2"),
  get_contributions_variables(species_single_comps$Resinicium_wl1) %>% mutate(comp = "Resinicium_wl1")
)

for_heatmap2 <- bind_rows(
  get_contributions_exp_variables(species_single_comps$Hypholoma_mill) %>% mutate(comp = "Hypholoma_mill"),
  get_contributions_exp_variables(species_single_comps$Hypholoma_nema) %>% mutate(comp = "Hypholoma_nema"),
  get_contributions_exp_variables(species_single_comps$Hypholoma_wl2) %>% mutate(comp = "Hypholoma_wl2"),
  
  get_contributions_exp_variables(species_single_comps$Phanerochaete_mill) %>% mutate(comp = "Phanerochaete_wl1"),
  # get_contributions_exp_variables(species_single_comps$Phanerochaete_wl2) %>% mutate(comp = "Phanerochaete_wl2"),
  # get_contributions_exp_variables(species_single_comps$Phanerochaete_nema) %>% mutate(comp = "Phanerochaete_nema"),
  
  get_contributions_exp_variables(species_single_comps$Resinicium_nema) %>% mutate(comp = "Resinicium_nema"),
  get_contributions_exp_variables(species_single_comps$Resinicium_wl2) %>% mutate(comp = "Resinicium_wl2"),
  get_contributions_exp_variables(species_single_comps$Resinicium_wl1) %>% mutate(comp = "Resinicium_wl1")
) %>% filter(grepl("Grazer", traits)) %>% group_by(comp) %>% summarise(RDA1 = sum(RDA1), RDA2 = sum(RDA2))


#for_heatmap2$comp <- factor(for_heatmap2$comp, levels = for_heatmap2$comp[order(for_heatmap2$RDA1)])

for_heatmap$traits <- factor(for_heatmap$traits, levels = scaled_variables)
for_heatmap$comp <- factor(for_heatmap$comp, levels = for_heatmap2$comp[order(for_heatmap2$RDA1)])
for_heatmap$comp2 <- factor(for_heatmap$comp, levels = for_heatmap2$comp[order(for_heatmap2$RDA2)])


cowplot::plot_grid(
  for_heatmap %>%
    ggplot() +
    aes(comp, fct_rev(traits), fill = RDA1) +
    geom_tile(color = "white", linewidth = 0.1) + 
    scale_fill_viridis_c(option = "plasma") +
    #scale_fill_gradient2(low = "red", mid = "white", high = "darkblue", midpoint = 0) +
    theme_minimal() +                          # Use a minimal theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis text
          axis.text.y = element_text(size = 8),
          axis.title = element_blank(),# Adjust y-axis text size
          panel.grid = element_blank()),
  
  
  for_heatmap %>%
    ggplot() +
    aes(comp2, fct_rev(traits), fill = RDA2) +
    geom_tile(color = "white", size = 0.1) +   # Add white borders to tiles
    scale_fill_viridis_c(option = "plasma") +  # Use 'plasma' color palette
    #scale_fill_gradient2(low = "red", mid = "white", high = "darkblue", midpoint = 0) +
    theme_minimal() +                          # Use a minimal theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis text
          axis.text.y = element_text(size = 8),
          axis.title = element_blank(),# Adjust y-axis text size
          panel.grid = element_blank())     
)


ggsave("heatmap_rda.png", plot = last_plot(), width = 14, height = 7, dpi = 300, units = "in")


all_data_scaled_a %>% 
  filter(!grepl("a_ctrl1|coll", Grazer)) %>% 
  group_by(Species, Grazer, day) %>%
  summarise(
    mean_angle = mean(Hyphal_angle),
    sd_angle = sd(Hyphal_angle)
  ) %>% 
  ggplot()+
  aes(x = day, y= mean_angle, colour = Grazer) +
  geom_point() +
    geom_path(aes(group = Grazer), linewidth = 2, arrow = arrow(), alpha = 0.7) +
    geom_ribbon(aes(ymin = mean_angle - sd_angle, ymax = mean_angle + sd_angle,
                    fill = Grazer, group = Grazer), alpha = 0.3) +
  scale_color_manual(
    values = c(
      "a_ctrl2" = "red",
      "mill" = "blue",
      "nema" = "green",
      "wl2" = "purple",
      "wl1" = "#D55E00"
      )
  ) +
  scale_fill_manual(
    values = c(
      "a_ctrl2" = "red",
      "mill" = "blue",
      "nema" = "green",
      "wl2" = "purple",
      "wl1" = "#D55E00"
    )
  ) +
  facet_wrap(.~Species, scales = "free")


###


all_data_scaled_a %>% 
  filter(grepl("Resi", Species)) %>% 
  filter(grepl("ctrl1|wl1", Grazer)) %>% 
  ggplot()+
  aes(x = day, y = Hyphal_angle) +
  geom_path(aes(color = Grazer, group = ID), linewidth = 2, arrow = arrow(), alpha = 0.7) +
  geom_text(aes(label = name_col)) +
  geom_point() 


contributions_exp_var <- bind_rows(
  get_contributions_exp_variables(species_single_comps$Hypholoma_mill) %>% mutate(comp = "Hypholoma_mill"),
  get_contributions_exp_variables(species_single_comps$Hypholoma_nema) %>% mutate(comp = "Hypholoma_nema"),
  get_contributions_exp_variables(species_single_comps$Hypholoma_wl2) %>% mutate(comp = "Hypholoma_wl2"),
  
  get_contributions_exp_variables(species_single_comps$Phanerochaete_wl1) %>% mutate(comp = "Phanerochaete_wl1"),
  get_contributions_exp_variables(species_single_comps$Phanerochaete_wl2) %>% mutate(comp = "Phanerochaete_wl2"),
  get_contributions_exp_variables(species_single_comps$Phanerochaete_nema) %>% mutate(comp = "Phanerochaete_nema"),
  
  get_contributions_exp_variables(species_single_comps$Resinicium_wl2) %>% mutate(comp = "Resinicium_wl2"),
  get_contributions_exp_variables(species_single_comps$Resinicium_wl1) %>% mutate(comp = "Resinicium_wl1")
)
  



### Supplementary figures

bind_rows(var_partiotion_single_comps, .id = "id") %>% 
  filter(!grepl("Total", Source)) %>%
  ggplot() +
  aes(x = Source, y = Proportion, fill = Source) + 
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(.~id, nrow = 4, ncol = 5) +
  theme(strip.text = element_text(size = 15))



