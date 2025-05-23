library(tidyverse)
library(readxl)
library(tidytext)
library(vegan)
library(ranger)

source("selected_load_data.R")

all_data_scaled_a %>% 
filter(day ==1) %>% 
  ggplot() +
  aes(PC1, PC2, colour = Species) +
  geom_point(size =2) +
  geom_text(aes(label = name_col))

#Experiment 1 (as October 2024, the only one)
exp1 <- grep("ctrl|coll|wl1|nema|wl1|wl2|mill", all_data_scaled_a$Grazer)

datos1 <- as_tibble(all_data_scaled_a[exp1, c("name_col", scaled_variables)]) %>% column_to_rownames(var = "name_col")

pca_exp1 <-
  rda(datos1[, scaled_variables],
      scale = TRUE, data = datos1)

all_interactions_var <- #This models assumes that each fungi has unique  and very idiosincratic development through time
  rda(all_data_scaled_a[exp1, scaled_variables]~ Species*Grazer*day + Condition(PC1 + PC2),
      scale = TRUE, data = all_data_scaled_a[exp1, ])

grazer_day_var <- # to me this makes the most sense. It assumes that each species has distinct starting point, but through time it has a common development. However, grazer can affect that development
  rda(all_data_scaled_a[exp1, scaled_variables]~ Species + Grazer + day + Condition(PC1 + PC2),
      scale = TRUE, data = all_data_scaled_a[exp1, ]);

# Based on this results, it looks like each species has its own trajectory. 
# This should be evident by comparing the trajectories of the controls across
# all species.

selected_model <-
  ordistep(grazer_day_var, scope = formula(all_interactions_var),
           permutations = 9999, direction = "forward")

selected_model_exp1 <- selected_model


#Getting summary stats
anova_exp1 <- 
  as.data.frame(
    anova.cca(selected_model_exp1, by = "term", permutations = 9999))

anova_exp1$R2 <- anova_exp1$Variance/sum(anova_exp1$Variance)

# The selected model (with all the interactions) explains:

summary(selected_model_exp1)[["cont"]][["importance"]]

######
# Measurment of effect sizes

# The code below with adonis2 does the same as anova.cca above but it less efficient (it takes longer to run it). The results are identical.
# The reason why I am repeating is because the custom made function "adonis_OmegaSq" to calculate effect sizes below
# only accepts adonis2 outputs.

# On October 2024, I had to stop using adonis2 approach because my new rda includes a Condition to take into account initial differences in 
# growth conditions. This "test" object I am making below is to confirm that the way R2 is measured in the adonis_OmegaSq function is simply
# a sum of the total variance and the variance explained by each term is divided the total sum. That is precisely the case so there is not need
# to use this approach really as the Omega sq was super correlated with R2 and R2 is easier to explain!
# 
# test <-  # to me this makes the most sense. It assumes that each species has distinct starting point, but through time it has a common development. However, grazer can affect that development
#   rda(all_data_scaled_a[exp1, scaled_variables]~ Species + Grazer,
#       scale = TRUE, data = all_data_scaled_a[exp1, ])
# 
# test <- as.data.frame(anova.cca(test, by = "term", permutations = 9999))
# 
# selected_model_exp1_ad <-
#   adonis2(scale(all_data_scaled_a[exp1, scaled_variables])~ Species + Grazer, # Based on the results of Ordistep
#           data = all_data_scaled_a[exp1,],
#           permutations = 9999,
#           method = "euclidean",
#           by = "terms",
#           na.rm = T) # It took 1 hour to run this!
# test$Variance/sum(test$Variance)
# 

#####

resultados_1 <-
 anova_exp1 %>% rownames_to_column(var = "Source")

resultados_1 <- resultados_1[, c(1:4, 6, 5)]

resultados_1 %>% 
  filter(grepl("Grazer",Source)) %>% 
  summarise(Source = "overall_grazer_effect", R2 = sum(R2))


write.table(rapply(resultados_1,
                   classes = "numeric", how = "replace",
                   round,digits = 4),
            "effect_sizes_full_model_exp1.txt", sep =";",
            row.names = F)



### Question 2. Which Grazer has the strongest impact on the similarity or dissmilarity of fungal colonies?


experiment_end <- all_data_scaled_a[, c("Species", "Grazer", "day", "name_col",
                                        "ID", "PC1", "PC2", scaled_variables)]
experiment_end <- experiment_end %>% filter(day == 8)

rda_grazer_function <- 
  function(dat, variables){rda(dat[,variables]~ Grazer + Condition (Species), scale = T, data = dat)}

grazer_single_comps <- 
  lapply(list(
    experiment_end %>%  filter(grepl("coll|a_ctrl", Grazer)), #
    experiment_end %>%  filter(grepl("wl1|a_ctrl", Grazer)), #
    experiment_end %>% filter(grepl("mill|a_ctrl", Grazer)), 
    experiment_end %>% filter(grepl("wl2|a_ctrl", Grazer)),
    experiment_end %>% filter(grepl("nema|a_ctrl", Grazer))),
    rda_grazer_function, variables = scaled_variables)

names(grazer_single_comps) <- c("collembola", "woodlice_1", "millipedes",
                                "woodlice_2", "nematodes")

anova_table <- function(model){
  as.data.frame(
    anova.cca(model, by = "term", permutations = 9999))
}

grazer_single_anova <- 
  lapply(grazer_single_comps, anova_table)

# Getting effect sizes

effect_size <- function(x){
  x$R2 <- x$Variance/sum(x$Variance)
  x}

grazer_effect_sizes <- 
  lapply(grazer_single_anova,
         effect_size)

table_grazer_effect_size <- 
  bind_rows(grazer_effect_sizes, .id = "id")

table_grazer_effect_size$source <- str_extract(rownames(table_grazer_effect_size), "[[:alpha:]]+")

#table_grazer_effect_size$source[seq(3, 20, by = 4)] <- "Species_x_Grazer"

rownames(table_grazer_effect_size) <- NULL
table_grazer_effect_size <- table_grazer_effect_size[,c(1,7,2:4,6,5)]
names(table_grazer_effect_size)[7] <- "p_value"

table_grazer_effect_size %>%
  filter(source == "Grazer")

table_grazer_effect_size %>% 
  filter(grepl("Grazer", source)) %>%
  filter(p_value < 0.06) %>%
  arrange(desc(R2)) %>% 
  ggplot() +
  aes(x = reorder(id, R2), y = R2) +
  geom_bar(stat = "identity") +
  #labs(title = "Deviation of treatment to control") +
  ylab(label = "Proportion of explained variance by grazer") +
  xlab(label = "Treatment") +
  theme(axis.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_blank()) +
  coord_flip()


write.table(rapply(table_grazer_effect_size[grep("Grazer|Species",
                                                 table_grazer_effect_size$source), ],
                   classes = "numeric", how = "replace",
                   round,digits = 4),
            "effect_sizes_per_grazer_sps_grazer.txt", sep =";",
            row.names = F)

write.table(rapply(table_grazer_effect_size,
                   classes = "numeric", how = "replace",
                   round,digits = 4),
            "effect_sizes_per_grazer_all.txt", sep =";",
            row.names = F)


# To complement the analysis, I will compare whether the dispersion between grazed and control (based on the means per treatments)

means_experiment_end <- 
experiment_end %>% 
  group_by(Species, Grazer) %>% 
  summarise_at(all_of(scaled_variables), mean)

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

# Function to calculate the PCA that takes into account the initial differences at the start of the experiment
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


# Function to extract the distances to the centroids as a dataframe
get_tables <- function(x){
  dat <- data.frame(treatments = x$group,
                    distances = x$distances)
  dat <- droplevels(dat)
}

# a) Calculate the PCA for each treatment

pca_grazers_wInitDiff <- 
lapply(list(
  experiment_end %>%  filter(grepl("coll|a_ctrl", Grazer)) %>% column_to_rownames(var = "name_col"), #
  experiment_end %>%  filter(grepl("wl1|a_ctrl", Grazer)) %>% column_to_rownames(var = "name_col"), #
  experiment_end %>% filter(grepl("mill|a_ctrl", Grazer)) %>% column_to_rownames(var = "name_col"), 
  experiment_end %>% filter(grepl("wl2|a_ctrl", Grazer)) %>% column_to_rownames(var = "name_col"),
  experiment_end %>% filter(grepl("nema|a_ctrl", Grazer)) %>% column_to_rownames(var = "name_col")),
  pca_function, variables = scaled_variables); names(pca_grazers_wInitDiff) <- c("collembola", "woodlice_1", "millipedes",
                                   "woodlice_2", "nematodes")


# b) Calculate the dispersion to each centroid and make a permutation test on that difference
s_grazer_dispersion_wInitDiff <- 
  lapply(pca_grazers_wInitDiff,
         dispersion_function, perm = T, group = T)

# c) Same as b, but without the test so I can extract the information to make the boxplot
s_plot_grazer_dispersion_wInitDiff <- 
  lapply(pca_grazers_wInitDiff,
         dispersion_function, perm = F, group = T)

lapply(s_plot_grazer_dispersion_wInitDiff, boxplot)

#For plotting

# I realize that the results of the stats above with the boxplot are qualitative the same as if I make one big
# pca (instead of separating each treatment). Doing one big PCA is easier to explai and to plot, so I decided to 
# go for that option




  

# What I decided is to show basic PC plot and the barplot, those ones show
# consistently the pattern. In the text I will explain that when testing these
# results using betadisper there was not statistical support for that difference
# this is in part because there was more "hidden" variation when one looks at 
# all the PC axes!
  
  

# Question 3. Which species is the most plastic and to what?

rda_grazer_sps_function <- 
  function(dat, variables){rda(dat[,variables]~ Grazer*day + Condition (PC1 + PC2), scale = T, data = dat)}# Orignal
  #function(dat, variables){rda(dat[,variables]~ day + Grazer:day , scale = T, data = dat)}# modified 28 August 2024
  #function(dat, variables){rda(dat[,variables]~ Grazer:day, scale = T, data = dat)}# modified on 25th August 2024


sps_experiment <- all_data_scaled_a[, c("Species", "Grazer", "day", "name_col", "ID",
                                        "PC1", "PC2", scaled_variables)]

######
# species_single_comps <- 
#   lapply(list(
#     #Hypholoma
#     sps_experiment %>% filter(grepl("coll|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("mill|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl1|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl2|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("nema|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),
#     
#     #Phallus
#     sps_experiment %>% filter(grepl("coll|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("mill|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl1|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl2|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("nema|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),
#     
#     #Phanerochaete
#     sps_experiment %>% filter(grepl("coll|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("mill|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl1|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl2|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("nema|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),
#     
#     #Resinicium
#     sps_experiment %>% filter(grepl("coll|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("mill|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl1|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("wl2|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col"),
#     sps_experiment %>% filter(grepl("nema|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col")),
#     
#     rda_grazer_sps_function, variables = scaled_variables)
# 
# names(species_single_comps) <- paste(
#   rep(c("Hypholoma", "Phallus", "Phanerochaete", "Resinicium"), each = 5),
#   rep(c("coll", "mill", "wl1", "wl2", "nema"), 4), sep = "_")

#####

species_single_comps <-
  lapply(list(
    #Hypholoma
    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Hypho", Species)) %>% column_to_rownames(var = "name_col"),

    #Phallus
    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Phallus", Species)) %>% column_to_rownames(var = "name_col"),

    #Phanerochaete
    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Phane", Species)) %>% column_to_rownames(var = "name_col"),

    #Resinicium

    sps_experiment %>% filter(grepl("mill|nema|wl2|wl1|coll|a_ctrl", Grazer)) %>% filter(grepl("Resi", Species)) %>% column_to_rownames(var = "name_col")),

    rda_grazer_sps_function, variables = scaled_variables)

names(species_single_comps) <- c("Hypholoma", "Phallus", "Phanerochaete", "Resinicium")


####

get_var_partition <- function(model){
  total_inertia <- model$tot.chi
  conditioned_inertia <- model$pCCA$tot.chi
  constrained_inertia <- model$CCA$tot.chi
  unconstrained_inertia <- model$CA$tot.chi
  
  # Calculate the proportions
  proportion_total <- 1
  proportion_conditioned <- conditioned_inertia / total_inertia
  proportion_constrained <- constrained_inertia / total_inertia
  proportion_unconstrained <- unconstrained_inertia / total_inertia
  
  # Create a data frame
  variance_partition <- data.frame(
    Inertia = c(conditioned_inertia, constrained_inertia, unconstrained_inertia,total_inertia),
    Proportion = c(proportion_conditioned, proportion_constrained, proportion_unconstrained, proportion_total),
    Source = c("Initial_diff", "Development_x_grazer", "Random", "Total")
  )
  variance_partition
}

var_partiotion_single_comps <-
lapply(species_single_comps, get_var_partition)

anova_table <- function(model){
  as.data.frame(
    anova.cca(model, by = "term", permutations = 9999))
}

species_single_anova <- 
  lapply(species_single_comps, anova_table)


# species_single_anova$source[seq(1, 40, by = 2)] <- paste(species_single_anova$source[seq(1, 40, by = 2)],
#                                                          "_x_day", sep = "")
# 
# species_single_anova$source[seq(2, 60, by = 3)] <- paste(species_single_anova$source[seq(2, 60, by = 3)],
#                                                          "_x_day", sep = "")

anova_table2 <- function(model){
  as.data.frame(
    anova.cca(model, by = "axis", permutations = 9999))
}

species_single_anova2 <- 
  lapply(species_single_comps, anova_table2)

# Getting effect sizes

species_single_anova <- 
  lapply(species_single_anova,
         effect_size)

species_single_anova <-
  bind_rows(species_single_anova, .id = "id")

species_single_anova$source <- 
  str_extract(rownames(species_single_anova), "[[:alpha:]]+")

species_single_anova$source[seq(3, length(species_single_anova$id), by = 4)] <- 
  paste(species_single_anova$source[seq(3, length(species_single_anova$id), by = 4)], "_x_day", sep = "")

rownames(species_single_anova) <- NULL

table_species_effect_size <- species_single_anova

table_species_effect_size <- table_species_effect_size[,c(1,7,2:4,6,5)]
names(table_species_effect_size)[7] <- "p_value"

saveRDS(table_species_effect_size, "comp_fungal_plasticity.RDS")

## Summing the effect of grazer and grazer x day

# table_species_effect_size %>% 
#   #filter(grepl("Grazer", source)) %>%
#   filter(p_value < 0.06 ) %>% 
#   group_by(id) %>% 
#   summarise(R2 = sum(R2)) %>% 
#   arrange(desc(R2)) %>% 
#   ggplot() +
#   aes(x = reorder(id, R2), y = R2) +
#   geom_bar(stat = "identity") +
#   #labs(title = "Deviation of treatment to control") +
#   ylab(label = "R2 of grazer and grazer:day") +
#   xlab(label = "Treatment") +
#   theme(axis.text = element_text(size = 30),
#         axis.title.x = element_text(size = 30)) +
#   coord_flip()

table_species_effect_size %>% 
  filter(grepl("Grazer", source)) %>%
  filter(p_value < 0.05 ) %>% 
  group_by(id) %>% 
  summarise(R2 = sum(R2)) %>% 
  arrange(desc(R2)) %>% 
  #filter(R2 > 0.1) %>% 
  ggplot() +
  aes(x = reorder(id, R2), y = R2) +
  geom_bar(stat = "identity") +
  #labs(title = "Deviation of treatment to control") +
  ylab(label = "Proportion of explained variation by Grazers") +
  theme(axis.text = element_text(size = 30),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 30)) +
  coord_flip()


ggsave("Figures//barplot_grazers_explained.png",
       plot = last_plot(), width = 14, height = 7, dpi = 300, units = "in")

# Possible supplementary figure

table_species_effect_size %>% 
  filter(grepl("_x", source)) %>%
  filter(p_value < 0.05) %>%
  arrange(desc(R2)) %>% 
  ggplot() +
  aes(x = reorder(id, R2), y = R2) +
  geom_bar(stat = "identity") +
  #labs(title = "Deviation of treatment to control") +
  ylab(label = "R2 of grazer:day only") +
  xlab(label = "Treatment") +
  theme(axis.text = element_text(size = 30),
        axis.title.x = element_text(size = 30)) +
  coord_flip()

write.table(rapply(table_species_effect_size %>% 
                     filter(p_value < 0.06) %>%
                     filter(grepl("Grazer_x", source)) %>% 
                     arrange(desc(R2)),
                   classes = "numeric", how = "replace",
                   round,digits = 4),
            "species_grazer_day_sig_effect_sizes.txt", sep =";",
            row.names = F)

write.table(rapply(table_species_effect_size,
                   classes = "numeric", how = "replace",
                   round,digits = 4),
            "species_grazer_day_all_effect_sizes.txt", sep =";",
            row.names = F)


