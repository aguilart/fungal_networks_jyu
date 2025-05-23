
library(tidyverse)
library(readxl)
library(tidytext)
library(vegan)
library(ranger)

rm(list = ls()); gc()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# Loading data

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

get_effect_size <- function(data){
  data <-  data %>%
      group_by(Grazer, day) %>%
      summarise(
        mean_PC1 = mean(PC1),
        var_PC1 = var(PC1),
        n = n(),
        mean_PC2 = mean(PC2),
        var_PC2 = var(PC2)
      ) %>% filter(day == 8)

cont <- which(data$Grazer=="z_ctrl")
coll <- which(data$Grazer=="coll")
mill <- which(data$Grazer=="mill")
nema <- which(data$Grazer=="nema")
wl1 <- which(data$Grazer=="wl1")
wl2 <- which(data$Grazer=="wl2")

data <- data.frame(Grazer = c("coll", "mill", "nema", "wl1", "wl2"),
                   t_test_PC1 = c(
(data$mean_PC1[cont] - data$mean_PC1[coll])/(sqrt((data$var_PC1[cont]/data$n[cont]) + (data$var_PC1[coll]/data$n[coll]))),
(data$mean_PC1[cont] - data$mean_PC1[mill])/(sqrt((data$var_PC1[cont]/data$n[cont]) + (data$var_PC1[mill]/data$n[mill]))),
(data$mean_PC1[cont] - data$mean_PC1[nema])/(sqrt((data$var_PC1[cont]/data$n[cont]) + (data$var_PC1[nema]/data$n[nema]))),
(data$mean_PC1[cont] - data$mean_PC1[wl1])/(sqrt((data$var_PC1[cont]/data$n[cont]) + (data$var_PC1[wl1]/data$n[wl1]))),
(data$mean_PC1[cont] - data$mean_PC1[wl2])/(sqrt((data$var_PC1[cont]/data$n[cont]) + (data$var_PC1[wl2]/data$n[wl2])))
                   ),
t_test_PC2 = c(
(data$mean_PC2[cont] - data$mean_PC2[coll])/(sqrt((data$var_PC2[cont]/data$n[cont]) + (data$var_PC2[coll]/data$n[coll]))),
(data$mean_PC2[cont] - data$mean_PC2[mill])/(sqrt((data$var_PC2[cont]/data$n[cont]) + (data$var_PC2[mill]/data$n[mill]))),
(data$mean_PC2[cont] - data$mean_PC2[nema])/(sqrt((data$var_PC2[cont]/data$n[cont]) + (data$var_PC2[nema]/data$n[nema]))),
(data$mean_PC2[cont] - data$mean_PC2[wl1])/(sqrt((data$var_PC2[cont]/data$n[cont]) + (data$var_PC2[wl1]/data$n[wl1]))),
(data$mean_PC2[cont] - data$mean_PC2[wl2])/(sqrt((data$var_PC2[cont]/data$n[cont]) + (data$var_PC2[wl2]/data$n[wl2])))
)

)}

prepare_data <- function(data) {
  # Summarize the original data for plotting
  plot_data <- data %>%
    group_by(Grazer, day) %>%
    summarise(
      mean_PC1 = mean(PC1),
      n = n(),
      sd_PC1 = sd(PC1)/sqrt(n),
      mean_PC2 = mean(PC2),
      sd_PC2 = sd(PC2)/sqrt(n)
    ) %>%
    pivot_longer(
      cols = starts_with("mean_") | starts_with("sd_"),
      names_to = c("stat", "variable"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = "stat",
      values_from = "value"
    ) %>%
    arrange(day)}

# Add a new column to determine which Grazer is highlighted
# plot_data <- plot_data %>%
#   mutate(highlight = ifelse(Grazer == highlight_grazer|Grazer == "z_ctrl", highlight_grazer, "Other"))


plot_pc_arrow_ribbons <- function(data, colores, var){
  # Plot with the specified Grazer highlighted
  data %>% filter(variable == var) %>% 
    ggplot() +
    aes(x = day, y = mean, group = Grazer) +
    geom_path(
      aes(color = Grazer),
      linewidth = 2, arrow = arrow(), alpha = 0.7
    ) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd, fill = Grazer),
      alpha = 0.3, outline.type = "both"
    ) +
    
    geom_line(
      aes(y = mean + sd, color = Grazer),
      linewidth = 1, alpha = 0.5
    ) +
    
    geom_line(
      aes(y = mean - sd, color = Grazer),
      linewidth = 1, alpha = 0.5
    ) +
    
    facet_grid(. ~ variable) +
    labs(title = unique(data$Species),
         x = "time",
         y = "PC") +
    scale_color_manual(
      values = colores
    ) +
    scale_fill_manual(
      values = colores
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",               # Legend at the bottom
      plot.title = element_text(size = 15))
}


###

scaled_variables <- c(#"Hyphal_main_width_median",
  
  "Hyphal_main_width_mean_of_MST", # To include only the main roads # 1
  "Hyphal_main_length_median_of_MST", # 2
  "Hyphal_sec_length_median", # 3
  "Hyphal_sec_width_mean",   # 4
  "Hyphal_angle", # 5
  "colony_Circularity", # 7
  
  
  # Measures of heterogeneity of the main roads
  "accessibility_skewness", # 8
  "length_skewness_main_of_MST", # 9
  "width_skewness_main_of_MST", # 10
  
 
  "Hyphal_density_main", # 6
  # Loopiness index
  "alpha_coeff", # 13
  "Mycelia_Vol_MST", # 14
  
  
  # Transport efficiency:
  "Reff_tip",  # 11
  #"Geff_MST",             
  "Geff_MST_no_c", # 12   
  
  
  
  # Robustness measures:
  
  "spatial_random_AUC_mean", # 16
  "spatial_random_AUC_cv" , # 17
  
  "AUC_Length_asc", # 18                   
  "AUC_Length_des", # 19
  "AUC_Resistance_asc", # 20               
  
  "AUC_Width_asc",  # 21                   
  "AUC_Width_des"#, # 22
)

###



location <- "C://Users//caaguila//Dropbox//Carlos Aguilar//MDF processing method comparison//"

# All methods deliver similar results, however after comaring them all visually and in pca the best ones were

# Resinicium SF0 M11 midgrey
# Hypholoma SF1 M11 midgrey
# Phanerochaete SF0 M11 midgrey
# Phallus SF0 M11 midgrey

# --------------- Loading data Hypholoma -----------------------------------

folder <- "Hf DD2//processed data SF1 M11 midgrey//results//"

Hf_T2_SF <- readRDS(paste(location, folder,"all_data_scaled.RDS", sep = ""))

Hf_T2_SF_mst <- readRDS(paste(location, folder,"all_data.RDS", sep = ""))
Hf_T2_SF_mst <- Hf_T2_SF_mst[which(Hf_T2_SF_mst$Network == "Resistance_MST"),
                             c("name_col", "summary_Geff",
                               "Geff_no_c",
                               "Root_eff",                         
                               "Reff_tip",
                               "accessibility_skewness",           
                               "length_skewness_tip",               "length_skewness_main",             
                               "volume_skewness_tip",               "volume_skewness_main",             
                               "width_skewness_tip",                "width_skewness_main",              
                               "Hyphal_tip_width_median",           "Hyphal_main_width_median",         
                               "Hyphal_main_width_mean", "Hyphal_main_length_median",
                               "Hyphal_tip_length_median")]
names(Hf_T2_SF_mst)[-1]  <- paste(names(Hf_T2_SF_mst)[-1], "_of_MST", sep = "")

Hf_T2_SF <- left_join(Hf_T2_SF, Hf_T2_SF_mst)

# adding new "colony" data also stored in the folder
temp = list.files(path = paste(location, folder, sep = ""), pattern = "*-results.xlsx")
temp <- paste(location, folder, temp,sep = "")

# For this I need to specify that the data is contained int an excel sheet called
# "Graph"
colony <- bind_rows(
  lapply(temp,function(x){read_excel(x,sheet = "Graph")} )); names(colony)[1] <- "name_col"

colony <- unique(colony[, c("name_col", "colony_Circularity")])
colony$name_col <- gsub("_shiny_", "_", colony$name_col)
Hf_T2_SF <- left_join(Hf_T2_SF, colony)

# Adding the new AUC measures
Carlos_random_rob <- readRDS(paste(location, folder,"Random_AUC.RDS", sep = ""))
names(Carlos_random_rob)[-1] <- paste("random_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_random_rob)[1]  <- "name_col"

Carlos_spatial_rob <- readRDS(paste(location, folder,"Spatial_AUC.RDS", sep = ""))
names(Carlos_spatial_rob)[-1] <- paste("spatial_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_spatial_rob)[1]  <- "name_col"

Carlos_ordered_rob <- readRDS(paste(location, folder,"ordered_AUC.RDS", sep = ""))
names(Carlos_ordered_rob)[1]  <- "name_col"

Carlos_rob <- 
  left_join(left_join(Carlos_ordered_rob, Carlos_spatial_rob),
            Carlos_random_rob)

#names(Carlos_rob)[-1] <- paste("CA_", names(Carlos_rob)[-1], sep = "")

Hf_T2_SF <- left_join(Hf_T2_SF,
                      Carlos_rob); rm(Carlos_rob, Carlos_ordered_rob, Carlos_spatial_rob, Carlos_random_rob)

Hf_T2_SF$set_up <- gsub("Hf DD2//|//results//|processed data ", "", folder)

# --------------- Loading data Resinicium -----------------------------------

folder <- "Rb//processed data SF0 M11 midgrey//results//"

Rb <- readRDS(paste(location, folder,"all_data_scaled.RDS", sep = ""))
Mst <- readRDS(paste(location, folder,"all_data.RDS", sep = ""))
Mst <- Mst[which(Mst$Network == "Resistance_MST"),
           c("name_col", "summary_Geff",
             "Geff_no_c",
             "Root_eff",                         
             "Reff_tip",
             "accessibility_skewness",           
             "length_skewness_tip",               "length_skewness_main",             
             "volume_skewness_tip",               "volume_skewness_main",             
             "width_skewness_tip",                "width_skewness_main",              
             "Hyphal_tip_width_median",           "Hyphal_main_width_median",         
             "Hyphal_main_width_mean", "Hyphal_main_length_median",
             "Hyphal_tip_length_median")]
names(Mst)[-1]  <- paste(names(Mst)[-1], "_of_MST", sep = "")
Rb <- left_join(Rb, Mst)

# adding new "colony" data also stored in the folder
temp = list.files(path = paste(location, folder, sep = ""), pattern = "*-results.xlsx")
temp <- paste(location, folder, temp,sep = "")

# For this I need to specify that the data is contained int an excel sheet called
# "Graph"
colony <- bind_rows(
  lapply(temp,function(x){read_excel(x,sheet = "Graph")} )); names(colony)[1] <- "name_col"

colony <- unique(colony[, c("name_col", "colony_Circularity")])
Rb <- left_join(Rb, colony)

# Adding Robustness

Carlos_random_rob <- readRDS(paste(location, folder,"Random_AUC.RDS", sep = ""))
names(Carlos_random_rob)[-1] <- paste("random_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_random_rob)[1]  <- "name_col"

Carlos_spatial_rob <- readRDS(paste(location, folder,"Spatial_AUC.RDS", sep = ""))
names(Carlos_spatial_rob)[-1] <- paste("spatial_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_spatial_rob)[1]  <- "name_col"

Carlos_ordered_rob <- readRDS(paste(location, folder,"ordered_AUC.RDS", sep = ""))
names(Carlos_ordered_rob)[1]  <- "name_col"

Carlos_rob <- 
  left_join(left_join(Carlos_ordered_rob, Carlos_spatial_rob),
            Carlos_random_rob)

#names(Carlos_rob)[-1] <- paste("CA_", names(Carlos_rob)[-1], sep = "")

Rb <- left_join(Rb, Carlos_rob); rm(Carlos_rob, Carlos_ordered_rob, Carlos_spatial_rob, Carlos_random_rob)

Rb$set_up <- gsub("Rb//|//results//|processed data ", "", folder)

# --------------- Loading data Phanaerochaete -----------------------------------

# Specifiying folder where data is located

folder <- "Pv//processed data SF0 M11 midgrey//results//"

# Loading the first type of data
Pv <- readRDS(paste(location, folder,"all_data_scaled.RDS", sep = ""))
Mst <- readRDS(paste(location, folder,"all_data.RDS", sep = ""))
Mst <- Mst[which(Mst$Network == "Resistance_MST"),
           c("name_col", "summary_Geff",
             "Geff_no_c",
             "Root_eff",                         
             "Reff_tip",
             "accessibility_skewness",           
             "length_skewness_tip",               "length_skewness_main",             
             "volume_skewness_tip",               "volume_skewness_main",             
             "width_skewness_tip",                "width_skewness_main",              
             "Hyphal_tip_width_median",           "Hyphal_main_width_median",         
             "Hyphal_main_width_mean", "Hyphal_main_length_median",
             "Hyphal_tip_length_median")]
names(Mst)[-1]  <- paste(names(Mst)[-1], "_of_MST", sep = "")
Pv <- left_join(Pv, Mst)

# adding new "colony" data also stored in the folder
temp = list.files(path = paste(location, folder, sep = ""), pattern = "*-results.xlsx")
temp <- paste(location, folder, temp,sep = "")

# For this I need to specify that the data is contained int an excel sheet called
# "Graph"
colony <- bind_rows(
  lapply(temp,function(x){read_excel(x,sheet = "Graph")} )); names(colony)[1] <- "name_col"

colony <- unique(colony[, c("name_col", "colony_Circularity")])
colony$name_col <- gsub("_shiny_", "_", colony$name_col)
Pv <- left_join(Pv, colony)

# Robustness

Carlos_random_rob <- readRDS(paste(location, folder,"Random_AUC.RDS", sep = ""))
names(Carlos_random_rob)[-1] <- paste("random_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_random_rob)[1]  <- "name_col"

Carlos_spatial_rob <- readRDS(paste(location, folder,"Spatial_AUC.RDS", sep = ""))
names(Carlos_spatial_rob)[-1] <- paste("spatial_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_spatial_rob)[1]  <- "name_col"

Carlos_ordered_rob <- readRDS(paste(location, folder,"ordered_AUC.RDS", sep = ""))
names(Carlos_ordered_rob)[1]  <- "name_col"

Carlos_rob <- 
  left_join(left_join(Carlos_ordered_rob, Carlos_spatial_rob),
            Carlos_random_rob)

#names(Carlos_rob)[-1] <- paste("CA_", names(Carlos_rob)[-1], sep = "")

Pv <- left_join(Pv, Carlos_rob); rm(Carlos_rob, Carlos_ordered_rob, Carlos_spatial_rob, Carlos_random_rob)

Pv$set_up <- gsub("Pv//|//results//|processed data ", "", folder)

# --------------- Loading data Phallus -----------------------------------

folder <- "Pi//processed data SF0 M11 midgrey//results//"

Pi <- readRDS(paste(location, folder,"all_data_scaled.RDS", sep = ""))
Mst <- readRDS(paste(location, folder,"all_data.RDS", sep = ""))
Mst <- Mst[which(Mst$Network == "Resistance_MST"),
           c("name_col", "summary_Geff",
             "Geff_no_c",
             "Root_eff",                         
             "Reff_tip",
             "accessibility_skewness",           
             "length_skewness_tip",               "length_skewness_main",             
             "volume_skewness_tip",               "volume_skewness_main",             
             "width_skewness_tip",                "width_skewness_main",              
             "Hyphal_tip_width_median",           "Hyphal_main_width_median",         
             "Hyphal_main_width_mean", "Hyphal_main_length_median",
             "Hyphal_tip_length_median")]
names(Mst)[-1]  <- paste(names(Mst)[-1], "_of_MST", sep = "")
Pi <- left_join(Pi, Mst)

# adding new "colony" data also stored in the folder
temp = list.files(path = paste(location, folder, sep = ""), pattern = "*-results.xlsx")
temp <- paste(location, folder, temp,sep = "")

# For this I need to specify that the data is contained int an excel sheet called
# "Graph"
colony <- bind_rows(
  lapply(temp,function(x){read_excel(x,sheet = "Graph")} )); names(colony)[1] <- "name_col"

colony <- unique(colony[, c("name_col", "colony_Circularity")])
Pi <- left_join(Pi, colony)

# Robustness
Carlos_random_rob <- readRDS(paste(location, folder,"Random_AUC.RDS", sep = ""))
names(Carlos_random_rob)[-1] <- paste("random_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_random_rob)[1]  <- "name_col"

Carlos_spatial_rob <- readRDS(paste(location, folder,"Spatial_AUC.RDS", sep = ""))
names(Carlos_spatial_rob)[-1] <- paste("spatial_", names(Carlos_random_rob)[-1], sep = "")
names(Carlos_spatial_rob)[1]  <- "name_col"

Carlos_ordered_rob <- readRDS(paste(location, folder,"ordered_AUC.RDS", sep = ""))
names(Carlos_ordered_rob)[1]  <- "name_col"

Carlos_rob <- 
  left_join(left_join(Carlos_ordered_rob, Carlos_spatial_rob),
            Carlos_random_rob)

#names(Carlos_rob)[-1] <- paste("CA_", names(Carlos_rob)[-1], sep = "")

Pi <- left_join(Pi, Carlos_rob); rm(Carlos_rob, Carlos_ordered_rob, Carlos_spatial_rob, Carlos_random_rob)

Pi$set_up <- gsub("Pi//|//results//|processed data ", "", folder)


# ------------------------------------------------------------------------------------------

# on top the densities and the widht, angle and AUC
library(vegan)
library(tidyverse)

all_data_scaled_c <- bind_rows(Hf_T2_SF , Pv, Rb, Pi)

all_data_scaled_c$Grazer <- NA
all_data_scaled_c$Grazer[grep("coll", all_data_scaled_c$name_col)] <- "coll"
all_data_scaled_c$Grazer[grep("ctrl1", all_data_scaled_c$name_col)] <- "ctrl1"
all_data_scaled_c$Grazer[grep("ctrl2", all_data_scaled_c$name_col)] <- "ctrl2"
all_data_scaled_c$Grazer[grep("mill", all_data_scaled_c$name_col)] <- "mill"
all_data_scaled_c$Grazer[grep("mite", all_data_scaled_c$name_col)] <- "mite"
all_data_scaled_c$Grazer[grep("nema", all_data_scaled_c$name_col)] <- "nema"
all_data_scaled_c$Grazer[grep("wl1", all_data_scaled_c$name_col)] <- "wl1"
all_data_scaled_c$Grazer[grep("wl2", all_data_scaled_c$name_col)] <- "wl2"

# There is grazer called "enchy" I remember we were going to exclude for several reasons
# I cannot remember them exactly. It only appears in Pi

all_data_scaled_c <-
  all_data_scaled_c[-which(is.na(all_data_scaled_c$Grazer)), ]

#all_data_scaled_c$name_col <- gsub(" ", "_" ,all_data_scaled_c$name_col)

all_data_scaled_c$day <- sapply(
  str_extract_all(all_data_scaled_c$name_col, "_D\\d+_|_d\\d+_"),
  function(x){x[1]})

all_data_scaled_c$day <- tolower(all_data_scaled_c$day)
all_data_scaled_c$day <- gsub("_", "", all_data_scaled_c$day)
all_data_scaled_c$day <- as.numeric(gsub("d", "", all_data_scaled_c$day))

all_data_scaled_c$ID <- sapply(
  strsplit(all_data_scaled_c$name_col, "_"),
  function(x){x[length(x)]})

all_data_scaled_c$ID <- paste(all_data_scaled_c$Species, all_data_scaled_c$Grazer,
                              all_data_scaled_c$ID, sep = "_")

all_data_scaled_c$ID <- gsub("Hypholoma fasiculare", "Hf",
                             all_data_scaled_c$ID)

all_data_scaled_c$ID <- gsub("Phallus impudicus", "Pi",
                             all_data_scaled_c$ID)

all_data_scaled_c$ID <- gsub("Phanaerochaete ventulina", "Pv",
                             all_data_scaled_c$ID)

all_data_scaled_c$ID <- gsub("Resinicium bicolor", "Rb",
                             all_data_scaled_c$ID)


all_data_scaled_c$Root_eff_a_scaled <- all_data_scaled_c$Root_eff/all_data_scaled_c$Mycelial_area
all_data_scaled_c$Reff_tip_a_scaled <- all_data_scaled_c$Reff_tip/all_data_scaled_c$Mycelial_area

#silence temporaly
#all_data_scaled_c$Geff_MST_no_c <- all_data_scaled_c$Geff_no_c/all_data_scaled_c$Geff_no_c_of_MST

# m <- which(all_data_scaled_c$day == 8|all_data_scaled_c$day == 1|all_data_scaled_c$day == 2|
#            all_data_scaled_c$day == 4|all_data_scaled_c$day == 16)
#all_data_scaled_a <- all_data_scaled_c[m,]

all_data_scaled_a <- all_data_scaled_c

all_data_scaled_a <- all_data_scaled_a[-grep("mite",
                                             all_data_scaled_a$Grazer), ]

table(all_data_scaled_a$Grazer,all_data_scaled_a$day, all_data_scaled_a$Species)


#Saving this file just to make it the loading faster
saveRDS(all_data_scaled_a, "current_all_data_scaled_a.RDS")
# Three replicates of Phanerochaete in Collembola are clear outliers and removed: Pv_coll_d1_3|Pv_col_d1_4|Pv_col_d1_5

#Clear outliers
all_data_scaled_a <- all_data_scaled_a[-grep("Pv_coll_d1_3|Pv_coll_d1_4|Pv_coll_d1_5", all_data_scaled_a$name_col), ]


start_exp <- all_data_scaled_a %>% 
  filter(day == 1)

first <- unique(start_exp$ID[which(start_exp$ID%in%all_data_scaled_a$ID)])

missing_first <-unique(
  all_data_scaled_a$ID[-which(all_data_scaled_a$ID%in%first)])

start_exp_2 <- all_data_scaled_a %>% 
  filter(ID%in%missing_first) %>% 
  filter(day == 2)

start_exp_all <- bind_rows(start_exp, start_exp_2)
rownames(start_exp_all) <- start_exp_all$ID

# initial_diff <- rda(start_exp_all[,scaled_variables], scale = T,
#                     data = start_exp_all)
# m <- sum(
#   length(initial_diff$CCA$eig),
#   length(initial_diff$CA$eig))
# test3 <- as.data.frame(scores(initial_diff,
#                               display = "sites",
#                               scaling = "sites",
#                               choices=c(1:m)))
# test3 <- rownames_to_column(test3) %>% rename(ID = rowname)

first_pca_function <- function(dat, variables){
  
  rda <- rda(dat[,variables], scale = T, data = dat)
  m <-
    sum(
      length(rda$CCA$eig),
      length(rda$CA$eig))
  
  test3 <- as.data.frame(scores(rda,
                                display = "sites",
                                scaling = "sites",
                                choices=c(1:m)))
  test3
}

initial_diff_list <- 
  lapply(list(
    #Hypholoma
    start_exp_all %>% filter(grepl("Hf", Species)),
    
    #Phallus
    start_exp_all %>% filter(grepl("Pi", Species)) ,
    
    #Phanerochaete
    start_exp_all %>% filter(grepl("Pv", Species)) ,
    
    #Resinicium
    start_exp_all %>% filter(grepl("Rb", Species))
  ),  first_pca_function, variables = scaled_variables)


test3 <- bind_rows(initial_diff_list)
test3 <- rownames_to_column(test3) %>% rename(ID = rowname)

###

all_data_scaled_a <- left_join(all_data_scaled_a, test3[,c("ID", "PC1", "PC2")],
                               by = "ID")

all_data_scaled_a$Grazer <- gsub("ctrl", "a_ctrl", all_data_scaled_a$Grazer)
all_data_scaled_a$Grazer[grep("ctrl", all_data_scaled_a$Grazer)] <- "a_ctrl"

all_data_scaled_a$Species[grep("Hf", all_data_scaled_a$Species)] <- "Hypholoma fasiculare"
all_data_scaled_a$Species[grep("Pv", all_data_scaled_a$Species)] <- "Phanerochaete velutina"
all_data_scaled_a$Species[grep("Rb", all_data_scaled_a$Species)] <- "Resinicium bicolor"
all_data_scaled_a$Species[grep("Pi", all_data_scaled_a$Species)] <- "Phallus impudicus"
all_data_scaled_a <- all_data_scaled_a %>% filter(day < 9)


# Functions

# Function to calculate differences in the dispersion of each grazer treatment to its centroid
dispersion_function <- function(traits, perm, group = F){
  
  if(group == T){
    traits <- traits %>% group_by(Species, Grazer) %>% 
      summarise_at(vars(PC1:PC21), mean)} else{traits}
  
  dis <- vegdist(traits[, grep("^PC", colnames(traits))],
                 na.rm = T, method = "euclidean")
  
  mod <- betadisper(dis, traits$Grazer)
  
  if(perm == T){ test <- permutest(mod, pairwise = TRUE, permutations = 9999)
  test}else{mod}
}


# Function to extract the distances to the centroids as a dataframe
get_tables <- function(x){
  dat <- data.frame(treatments = x$group,
                    distances = x$distances)
  dat <- droplevels(dat)
}