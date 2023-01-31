### Around the island nutrient regimes ###
### Script by Nyssa Silbiger ###
### Updated on 2023-01-30 #######


##### Load Libraries #############
library(tidyverse)
library(here)
library(biscale)
library(patchwork)
library(cowplot)



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
