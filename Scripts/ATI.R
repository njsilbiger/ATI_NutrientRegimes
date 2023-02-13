### Around the island nutrient regimes ###
### Script by Nyssa Silbiger ###
### Updated on 2023-01-30 #######


##### Load Libraries #############
library(tidyverse)
library(here)
library(biscale)
library(patchwork)
library(cowplot)
library(PNWColors)
library(hrbrthemes)
library(circlize)
library(networkD3)
library(ggsankey)
library(ggalluvial)
library(FactoMineR)
library(ggfortify)
library(cluster)

#### read in the data ###############
data<-read_csv(here("Data","AllNutrientData_clusters.csv"))


# filter to just the sites and times with water column nutrient data ####

data_filtered <-data %>%
  drop_na(Nitrite_plus_Nitrate) %>%
  filter(Year != 2020)  %>% # drop the bad 2020 data
  filter(!is.na(Turbinaria_groups)) # remove the data that have no turbinaria groups

### scale the data by year ####
data_scaled<-data_filtered %>%
  group_by(Year) %>%
  mutate_at(vars(c("Phosphate","Silicate","Nitrite_plus_Nitrate","Ammonia",
                   "Ultra.Violet.Humic.like","Tyrosine.like","Visible.Humic.like",
                   "Marine.Humic.like","Tryptophan.like","Lignin.like","M.C",
                   "HIX","BIX","FI")), function(x){log(x+0.1)})%>%
  mutate_at(vars(c("Phosphate","Silicate","Nitrite_plus_Nitrate","Ammonia",
                   "Ultra.Violet.Humic.like","Tyrosine.like","Visible.Humic.like",
                   "Marine.Humic.like","Tryptophan.like","Lignin.like","M.C",
                   "HIX","BIX","FI")), list(scale = scale) )%>%
  ungroup()


### run a PCA on the nutrient data with centered data by year #####

## Just the nutrients ###

nutrientdata<-data_scaled %>%
  select(Phosphate_scale:Ammonia_scale)

pca_nuts<-princomp(nutrientdata)

data_scaled_pca<-bind_cols(data_scaled, data.frame(pca_nuts$scores[,1:2]))

## Choose bivariate color palette and generate legend
bivColPal <- "GrPink"
legendRas <- bi_legend(pal = bivColPal,
                       dim = 3,
                       xlab = "Higher % N",
                       ylab = "Higher variance ",
                       size = 12)
# add the color pallete
data_scaled_pca <- biscale::bi_class(data_scaled_pca, 
                                 x = Turbinaria_N, 
                                 y = Turbinaria_V,
                                 style = "quantile", 
                                 dim = 3)

p1<-ggplot(data_scaled_pca)+
  geom_point(aes(x = Comp.1, y = Comp.2, color = bi_class),
             show.legend = FALSE)+
  facet_wrap(~Year) +
  theme_bw()+
  bi_scale_color(pal = bivColPal, dim = 3) +
  ggnewscale::new_scale("color") +
  theme(line = element_blank())


# make biplot
loadingdata<-data.frame(pca_nuts$loadings[,1:2])

#p1+legendRas

finalPlot <- cowplot::ggdraw() +
  cowplot::draw_plot(p1, 0, 0, 1, 1)+
  cowplot::draw_plot(legendRas, 0.025, 0.65, 0.275, 0.3) 
  
biplot <- ggplot(loadingdata)+
  coord_equal() + 
  geom_text(aes(x=Comp.1, y=Comp.2, label=rownames(loadingdata)), size = 5, vjust=1, color="red")+
  geom_segment(aes(x=0, y=0, xend=Comp.1, yend=Comp.2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")+
  theme_bw()+
  theme(line = element_blank())

finalPlot|biplot
ggsave(here("Outputs","nutrientpca.png"), width = 8, height = 6)

### Same with fDOM included


nutrientfdom_data<-data_scaled %>%
  select(Phosphate_scale:Lignin.like_scale) %>%
  drop_na(Lignin.like_scale)

pca_nuts_fdom<-princomp(nutrientfdom_data)

data_scaled_pca2<-bind_cols(data_scaled%>%
                              drop_na(Lignin.like_scale), data.frame(pca_nuts_fdom$scores[,1:2]))

# add the color pallete
data_scaled_pca2 <- biscale::bi_class(data_scaled_pca2, 
                                     x = Turbinaria_N, 
                                     y = Turbinaria_V,
                                     style = "quantile", 
                                     dim = 3)

p2<-ggplot(data_scaled_pca2)+
  geom_point(aes(x = Comp.1, y = Comp.2, color = bi_class),
             show.legend = FALSE)+
  facet_wrap(~Year) +
  theme_bw()+
  bi_scale_color(pal = bivColPal, dim = 3) +
  ggnewscale::new_scale("color") +
  theme(line = element_blank())

# make biplot
loadingdata2<-data.frame(pca_nuts_fdom$loadings[,1:2])

biplot2 <- ggplot(loadingdata2)+
  coord_equal() + 
  geom_text(aes(x=Comp.1, y=Comp.2, label=rownames(loadingdata2)), size = 5, vjust=1, color="red")+
  geom_segment(aes(x=0, y=0, xend=Comp.1, yend=Comp.2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")+
  theme_bw()+
  theme(line = element_blank())

# bring it with the legend
finalPlot2 <- cowplot::ggdraw() +
  cowplot::draw_plot(p2, 0, 0, 1, 1)+
  cowplot::draw_plot(legendRas, 0.725, 0.1, 0.275, 0.3) 

fdom<-finalPlot2|biplot2

ggsave(filename = here("Outputs","nutrientfdompca.png"), plot = fdom, width = 12, height = 6)

# Nury's clusters
set.seed(123)
fviz_nbclust(nutrientfdom_data, pam, method="silhouette")+theme_classic()

set.seed(123)
pam.res <- pam(nutrientfdom_data, 3) #using 5 even tho the silhoutte chose 2
print(pam.res)

## add in Nury's clusters
data_scaled_pca2  <- data_scaled_pca2 %>%
  mutate(Cluster_Nury = factor(pam.res$cluster))

# add Nury's clusters with the mega dataset
data<-data_scaled_pca2 %>%
  select(Site, Year, Cluster_Nury) %>%
  right_join(data)

# make a PCA with hers

ggplot(data_scaled_pca2)+
  geom_point(aes(x = Comp.1, y = Comp.2, color = Cluster_Nury),
             show.legend = FALSE)+
  facet_wrap(~Year) +
  theme_bw()+
  theme(line = element_blank())
ggsave(here("Outputs","Nury_clusters.png"), width = 6, height = 4)
### Calculate the euclidean distances between the years for the pca

distance<- function(x1, x2, y1, y2){
  sqrt((x2-x1)^2+(y2-y1)^2)
  }

Y21<-data_scaled_pca2 %>%
  select(Site, Year, Turbinaria_groups, Comp.1, Comp.2, Nutrient_Clusters) %>%
  filter(Year == "2021") %>%
  select(-Year, Comp.1_21 = Comp.1, Comp.2_21 = Comp.2, Nutrient_Clusters)

Y22<-data_scaled_pca2 %>%
  select(Site, Year, Turbinaria_groups, Comp.1, Comp.2, Nutrient_Clusters) %>%
  filter(Year == "2022") %>%
select(-Year, Comp.1_22 = Comp.1, Comp.2_22 = Comp.2, Nutrient_Clusters)

distances<-full_join(Y21, Y22)  %>%
  mutate(Eu_distance = distance(Comp.1_21, Comp.1_22, Comp.2_21, Comp.2_22))

distances %>%
  ggplot(aes(x = Turbinaria_groups, y = Eu_distance))+
  geom_boxplot()

# color palette
pal <- pnw_palette("Sailboat",3, type = "discrete")

distances %>%
  drop_na()%>%
  ggplot(aes(x = Nutrient_Clusters, y = Eu_distance, fill = Nutrient_Clusters))+
  # ggdist::stat_gradientinterval(
  #   linewidth = .3, color = "black"
  # )+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  labs(x = "",
       y = "Euclidean distance between 2021 and 2022",
       fill = "Nutrient Clusters")+
  scale_fill_manual(values = pal)+
  theme_bw()
ggsave(here("Outputs","Eu_distances.png"), width = 6, height = 6)

# try with Craigs PCA
Y21a<-data %>%
  select(Site, Year, Turbinaria_groups, Nutrient_PC1, Nutrient_PC2, Nutrient_Clusters, Cluster_Nury) %>%
  filter(Year == "2021") %>%
  select(-Year, Comp.1_21 = Nutrient_PC1, Comp.2_21 = Nutrient_PC2, Nutrient_Clusters, Cluster_Nury)

Y22a<-data %>%
  select(Site, Year, Turbinaria_groups, Nutrient_PC1, Nutrient_PC2, Nutrient_Clusters,Cluster_Nury) %>%
  filter(Year == "2022") %>%
  select(-Year, Comp.1_22 = Nutrient_PC1, Comp.2_22 = Nutrient_PC2, Nutrient_Clusters2 = Nutrient_Clusters,Cluster_Nury2 = Cluster_Nury)

distancesa<-full_join(Y21a, Y22a)  %>%
  mutate(Eu_distance = distance(Comp.1_21, Comp.1_22, Comp.2_21, Comp.2_22))

# figure out the percent of points that changed clusters
distancesa %>%
  mutate(same = ifelse(Nutrient_Clusters == Nutrient_Clusters2, "yes","no")) %>%
  drop_na(Nutrient_Clusters,Nutrient_Clusters2)%>%
  group_by(same) %>%
  count()

# 60% stayed in the same cluster

distancesa %>%
  drop_na()%>%
  ggplot(aes(x = Nutrient_Clusters, y = Eu_distance))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(color = Nutrient_Clusters2))+
  labs(x = "Nutrient Clusters 2021",
       y = "Euclidean distance between 2021 and 2022",
       color = "Nutrient Clusters 2022")+
  scale_color_manual(values = pal)+
  theme_bw()

ggsave(here("Outputs","Eu_distances_craig.png"), width = 6, height = 6)

## make a sankey plot
# make the df with counts for each group
sankeydf<-distancesa %>%
  drop_na()%>%
  mutate(same = ifelse(Nutrient_Clusters == Nutrient_Clusters2, "yes","no"))%>%
  group_by(Nutrient_Clusters, Nutrient_Clusters2, same) %>%
  count(.drop = FALSE) %>%
  ungroup()%>% # calculate the number within each group
  select(source = Nutrient_Clusters,target = Nutrient_Clusters2,  value = n, same)%>%
  bind_rows(data.frame(source = c("High Nutrient"), target = c("High fDOM"), value = c(0), same = c("no")))

# make the plot
ggplot(data = sankeydf,
       aes(axis1 = source, axis2 = target, 
           y = value)) +
  scale_x_discrete(limits = c("2021", "2022"), expand = c(.2, .05)) +
 # xlab("Demographic") +
  geom_alluvium(aes(fill = source)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+
  scale_fill_manual(values = pal)+
  labs(y = "",
      # fill = "Stayed the same"
       )+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')

ggsave(here("Outputs","flowdiagram.jpg"), width = 6, height = 6)

# make the plot
ggplot(data = sankeydf2,
       aes(axis1 = source, axis2 = target, 
           y = value)) +
  scale_x_discrete(limits = c("2021", "2022"), expand = c(.2, .05)) +
  # xlab("Demographic") +
  geom_alluvium(aes(fill = source)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+
  scale_fill_manual(values = pal)+
  labs(y = "",
       # fill = "Stayed the same"
  )+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')
## Same with Nury clusters
sankeydf2<-distancesa %>%
  drop_na()%>%
  mutate(same = ifelse(Cluster_Nury == Cluster_Nury2, "yes","no"))%>%
  group_by(Cluster_Nury, Cluster_Nury2, same) %>%
  count(.drop = FALSE) %>%
  ungroup()%>% # calculate the number within each group
  select(source = Cluster_Nury,target = Cluster_Nury2,  value = n, same)

ggsave(here("Outputs","Nuryclusterschange.jpg"), width = 6, height = 6)
## plots of craig's PC2 versus mean turb
data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Turbinaria_quantile_meanName, y = Nutrient_PC2, fill = Turbinaria_quantile_meanName))+
  geom_boxplot()+
  geom_jitter(width = 0.2, alpha = 0.5)+
  labs(x = "Turbinarea mean",
       y = "PC2",
       fill = "Turb Mean Clusters")+
  scale_fill_manual(values = pal)+
  theme_bw()

mod_N<-lm(Nutrient_PC2~Turbinaria_N, data = data %>% drop_na(Turbinaria_N,Nutrient_PC2))
anova(mod_N)
summary(mod_N)

mod_N_type2<-lmodel2(Nutrient_PC2~Turbinaria_N, data = data %>% drop_na(Turbinaria_N,Nutrient_PC2))

fits<-tidy(mod_N_type2)
slope <-fits$estimate[6]
intercept <-fits$estimate[5]

data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Percent_N, y = Nutrient_PC2))+
  geom_point(aes(color = factor(Year)))+
  geom_smooth(method = "lm")+
  geom_abline(slope = slope, intercept = intercept, color = "black")+
  labs(color = "Year",
      x = "Turbinatia %N Mean",
      y = "Nutrient PC2")+
  theme_bw()

ggsave(here("Outputs","Turb_pc2.png"), width = 6, height = 6)

mod_N<-lm(Nutrient_PC2~Turbinaria_N, data = data %>% drop_na(Turbinaria_N,Nutrient_PC2))
anova(mod_N)
summary(mod_N)

mod_N_type2<-lmodel2(Nutrient_PC2~Turbinaria_N, data = data %>% drop_na(Turbinaria_N,Nutrient_PC2))
summary(mod_N_type2)

mod_N2<-lm(Nutrient_PC2~Percent_N, data = data %>% drop_na(Turbinaria_N,Nutrient_PC2))
anova(mod_N2)
summary(mod_N2)

# Variance model
data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Turbinaria_V, y = Nutrient_PC2))+
  geom_point(aes(color = factor(Year)))+
  geom_smooth(method = "lm")+
  labs(color = "Year",
       x = "Turbinatia %N Variance",
       y = "Nutrient PC2")+
  theme_bw()
ggsave(here("Outputs","Turb_pc2_var.png"), width = 6, height = 6)

## Nury clusters with microbes
N1<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Cluster_Nury, y = Microbial_PCoA1, fill = Cluster_Nury))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  scale_fill_manual(values = pal)+
 # geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial PCO1",
    x = "Nutrient cluster Nury")+
  theme_bw()+
  theme(legend.position = 'none')

N2<-data %>%
  mutate(Nutrient_Clusters = factor(Nutrient_Clusters, levels = c("Low Nutrient","High fDOM","High Nutrient")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_Clusters, y = Microbial_PCoA1, fill = Nutrient_Clusters))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  scale_fill_manual(values = pal)+
  # geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial PCO1",
    x = "Nutrient cluster Craig")+
  theme_bw()+
  theme(legend.position = 'none')

N1+N2
ggsave(here("Outputs","clustercomparison.png"), width = 8, height = 4)
### with the microbe data
PM1<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC2, y = Microbial_PCoA1))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
       y = "Microbial PCO1",
       x = "Nutrient cluster PC2")+
  theme_bw()

MicroPC1mod<-lm(Microbial_PCoA1~Nutrient_PC2, data = data)
summary(MicroPC1mod)


PM2<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC2, y = Microbial_Species_Richness))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Species Richness",
    x = "Nutrient cluster PC2")+
  theme_bw()

MicroRichmod<-lm(Microbial_Species_Richness~Nutrient_PC2, data = data)
summary(MicroRichmod)

PM3<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC2, y = Microbial_Shannon_Diversity))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Shannon Diversity",
    x = "Nutrient cluster PC2")+
  theme_bw()

MicroShanmod<-lm(Microbial_Shannon_Diversity~Nutrient_PC2, data = data)
summary(MicroShanmod)


PM4<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC2, y = Microbial_Phylogenetic_Diversity))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Phylogenetic Diversity",
    x = "Nutrient cluster PC2")+
  theme_bw()

MicroPhylomod<-lm(Microbial_Phylogenetic_Diversity~Nutrient_PC2, data = data)
summary(MicroPhylomod)

PM5<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC2, y = Microbial_Evenness))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Evenness",
    x = "Nutrient cluster PC2")+
  theme_bw()

MicroEvenmod<-lm(Microbial_Evenness~Nutrient_PC2, data = data)
summary(MicroEvenmod)

PM6<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC2, y = Microbial_PCoA2))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial PCOA2",
    x = "Nutrient cluster PC2")+
  theme_bw()

MicroPCO2mod<-lm(Microbial_PCoA2~Nutrient_PC2, data = data)
summary(MicroPCO2mod)


(PM1 +PM2+PM4)/ (PM6+PM3+PM5)
ggsave(here("Outputs","Microbe_PC2.png"), width = 8, height = 6)

################ same with PC cluster 1
PM1a<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC1)%>%
  ggplot(aes(x = Nutrient_PC1, y = Microbial_PCoA1))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial PCO1",
    x = "Nutrient cluster PC1")+
  theme_bw()

PM2a<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC1, y = Microbial_Species_Richness))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Species Richness",
    x = "Nutrient cluster PC1")+
  theme_bw()

PM3a<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC1, y = Microbial_Shannon_Diversity))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Shannon Diversity",
    x = "Nutrient cluster PC1")+
  theme_bw()

PM4a<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC1, y = Microbial_Phylogenetic_Diversity))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Phylogenetic Diversity",
    x = "Nutrient cluster PC1")+
  theme_bw()

PM5a<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC2)%>%
  ggplot(aes(x = Nutrient_PC1, y = Microbial_Evenness))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial Evenness",
    x = "Nutrient cluster PC1")+
  theme_bw()

PM6a<-data %>%
  mutate(Turbinaria_quantile_meanName = factor(Turbinaria_quantile_meanName, levels = c("Low","Med","High")))%>%
  drop_na(Turbinaria_N,Nutrient_PC1)%>%
  ggplot(aes(x = Nutrient_PC1, y = Microbial_PCoA2))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(#color = "Year",
    y = "Microbial PCOA2",
    x = "Nutrient cluster PC1")+
  theme_bw()

(PM1a +PM2a+PM4a)/ (PM6a+PM3a+PM5a)
