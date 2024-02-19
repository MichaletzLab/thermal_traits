###...................................###
# Bison NN, Michaletz ST. Variation in leaf carbon economics, energy balance, and heat tolerance traits highlights differing time scales of adaptation and acclimation 
# Code S1: R code for reproducing analyses and figures. 
# Prepared by Nicole Bison (nicole.bison@ubc.ca), September 2023, Revised February 2024.
###...................................###

#Load packages
library(htmltools)
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(tidyverse)
library(V.PhyloMaker)
library(plantlist)
library(dplyr)
library(ggtree)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(ggjoy)
#library("devtools")
library(vctrs)
library(htmltools)
#library("FactoMineR")
library("factoextra")
#install.packages("hier.part")
library(hier.part)
library(moments)
library(caper)
library(rr2)
library(phylolm)
library(caper)

# Load dataset with taxonomic information & trait data
phy_traits_sp_raw <- read.csv("data_S1_revised_v3.csv")
phy_traits_sp <- phy_traits_sp_raw

########### Main text ###########
###### Figure 2: Trait distributions  ####
tcrit_ridge <- ggplot(phy_traits_sp, aes(x = Tcrit, y = clade, fill = clade, color = clade)) +
  scale_fill_manual(values = c("#40004B", "#762A83","#9970AB"))+
  scale_colour_manual(values = c("#40004B", "#762A83","#9970AB"))+
  xlab("Tcrit (°C)")+
  scale_x_continuous(trans = 'log10')+
  ylab("")+
  xlab(expression(T[crit]~("°C"))) +
  ylab("")+
  guides(fill = guide_legend(reverse=TRUE),color = guide_legend(reverse=TRUE))+
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'), 
        legend.box.background = element_rect(fill='transparent')) +
  theme(legend.background = element_rect(color = NA))+
  geom_density_ridges(alpha = 0.8)
tcrit_ridge

tau_ridge<- ggplot(phy_traits_sp, aes(x = tau, y = clade, fill = clade, color = clade)) +
  scale_colour_manual(values = c("#21484A","#336568","#488087"))+
  scale_fill_manual(values = c("#21484A","#336568","#488087"))+
  labs(x = tau~"(s)") +
  ylab("Density")+
  scale_x_continuous(trans = 'log10')+
  guides(fill = guide_legend(reverse=TRUE),color = guide_legend(reverse=TRUE))+
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'), 
        legend.box.background = element_rect(fill='transparent')) +
  theme(legend.background = element_rect(color = NA))+
  geom_density_ridges(alpha = 0.8)
tau_ridge

LMA_ridge<- ggplot(phy_traits_sp, aes(x = lma, y = clade, fill = clade, color = clade)) +
  scale_fill_manual(values = c("#6C2C15", "#8D4A21", "#AD692D"))+
  scale_colour_manual(values = c("#6C2C15", "#8D4A21", "#AD692D"))+
  scale_x_continuous(trans = 'log10')+
  xlab(bquote("LMA (kg "~m^-2 ~ ")"))+
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA), 
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  ylab("")+  
  guides(fill = guide_legend(reverse=TRUE),color = guide_legend(reverse=TRUE))+
  geom_density_ridges(alpha = 0.8)
LMA_ridge

LDMC_ridge<- ggplot(phy_traits_sp, aes(x = ldmc, y = clade, fill = clade, color = clade)) +
  scale_fill_manual(values = c("#6C2C15", "#8D4A21", "#AD692D"))+
  scale_colour_manual(values = c("#6C2C15", "#8D4A21", "#AD692D"))+
  scale_x_continuous(trans = 'log10')+
  xlab(bquote("LDMC (kg "~kg^-1 ~ ")"))+
  ylab("")+
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA), 
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  theme(legend.background = element_rect(color = NA))+
  ylab("")+  
  guides(fill = guide_legend(reverse=TRUE),color = guide_legend(reverse=TRUE))+
  geom_density_ridges(alpha = 0.8)
LDMC_ridge



ridge_plot_supp <- ggarrange(LMA_ridge, LDMC_ridge, tau_ridge,tcrit_ridge,
                             labels = c("a", "b","c", "d"), hjust = -4, vjust =1, 
                             ncol = 1, nrow = 4)
ridge_plot_supp
ggsave("./figures_revised/Figure2.pdf", plot = ridge_plot_supp, width = 17, height = 22, units = "cm", dpi = 300)


###### Table 2: Test for phylogenetic signal ######
#Build tree from list of species names
sp_list <- phy_traits_sp[c(3,4,2)]

#Log transform traits to address non-normality
phy_traits_sp$Tcrit <- log(phy_traits_sp$Tcrit)
phy_traits_sp$Tmax <- log(phy_traits_sp$Tmax)
phy_traits_sp$lma <- log(phy_traits_sp$lma)
phy_traits_sp$ldmc <- log(phy_traits_sp$ldmc)
phy_traits_sp$tau <- log(phy_traits_sp$tau)
phy_traits_sp$char_dimension <- log(phy_traits_sp$char_dimension)
phy_traits_sp$leaf_area <- log(phy_traits_sp$leaf_area)

tree <- phylo.maker(sp_list, tree = GBOTB.extended, output.tree = T, r = 1)
tree3 <- tree$scenario.3


# Also prep dataset for PGLS, so that we can further test whether controlling for seasonality in dataset impacts phylogenetic signal
tree_up <- tree3
tree_up$node.label <- NULL
comp.data<-comparative.data(tree_up,phy_traits_sp, "name")

#Tcrit
#Using phylosig, no accounting for sampling date
trait <- setNames(phy_traits_sp$Tcrit, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

#Using PGLS to account for when that species was sampling during the study period (mean temperature on sampling date)
model<-pgls(Tcrit  ~  mean_temp, data=comp.data, lambda = "ML")
summary(model)

#Tmax
#Using phylosig, no accounting for sampling date
trait <- setNames(phy_traits_sp$Tmax, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

#Using PGLS to account for when that species was sampling during the study period (mean temperature on sampling date)
model<-pgls(Tmax ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

#lma
trait <- setNames(phy_traits_sp$lma, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(lma ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

#ldmc
trait <- setNames(phy_traits_sp$lma, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(ldmc ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

#tau
trait <- setNames(phy_traits_sp$tau, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(tau ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

#area
trait <- setNames(phy_traits_sp$leaf_area, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(leaf_area ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

#characteristic dimension
trait <- setNames(phy_traits_sp$char_dimension, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(char_dimension ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

###### Figure 3: Plotting traits on phylogeny ######
#Calculate quantiles for each trait
#Not using log transformed values for plotting so we can see the actual range for measured values
#Binning  traits to address long tails of distribution in colour scheme
phy_traits_sp_plotting <- phy_traits_sp_raw

tcrit_quantiles <- quantile(phy_traits_sp_plotting$Tcrit, prob=c(0.05,.25,.50,.75,.95), type=1)
tau_quantiles <-quantile(phy_traits_sp_plotting$tau, prob=c(0.05,.25,.50,.75,.95), type=1)
lma_quantiles <- quantile(phy_traits_sp_plotting$lma, prob=c(0.05,.25,.50,.75,.95), type=1)
ldmc_quantiles <- quantile(phy_traits_sp_plotting$ldmc, prob=c(0.05,.25,.50,.75,.95), type=1)
abs_quantiles <- quantile(phy_traits_sp_plotting$absorptance, prob=c(0.05,.25,.50,.75,.95), type=1, na.rm = T)

#Create bins with the quantiles calculated above
df_d <- phy_traits_sp_plotting
df_d$name <- sub(" ", "_", df_d$name)

df1 <- as.data.frame(df_d$tau)
colnames(df1)[1] <- "tau"
tau_quantiles <- unname(tau_quantiles)
df1 <- df1 %>% mutate(tau_bin = cut(tau, breaks=c(min(phy_traits_sp_raw$tau) - 0.001,tau_quantiles[1],tau_quantiles[2], #-0.01 ensures min observation fits in bin
                                                  tau_quantiles[3],tau_quantiles[4],tau_quantiles[5] ,max(phy_traits_sp_raw$tau + 0.001))))  #+0.01 ensures max observation fits in bin
df1 <- as.data.frame(df1[,2])
colnames(df1)[1] <- "tau_bin"
rownames(df1) <- df_d$name


df2 <- as.data.frame(df_d$Tcrit)
colnames(df2)[1] <- "Tcrit"
tcrit_quantiles <- unname(tcrit_quantiles)
df2 <- df2 %>% mutate(tcrit_bin = cut(Tcrit, breaks=c(min(phy_traits_sp_raw$Tcrit)-0.001,tcrit_quantiles[1],tcrit_quantiles[2],
                                                      tcrit_quantiles[3],tcrit_quantiles[4],tcrit_quantiles[5] ,max(phy_traits_sp_raw$Tcrit)+0.001)))
df2 <- as.data.frame(df2[,2])
colnames(df2) <- "tcrit_bin"
rownames(df2) <- df_d$name


df3 <- as.data.frame(df_d$lma)
colnames(df3)[1] <- "lma"
lma_quantiles <- unname(lma_quantiles)
df3 <- df3 %>% mutate(lma_bin = cut(lma, breaks=c(min(phy_traits_sp_raw$lma)-0.001,lma_quantiles[1],lma_quantiles[2],
                                                  lma_quantiles[3],lma_quantiles[4],lma_quantiles[5] ,max(phy_traits_sp_raw$lma + 0.001))))
df3 <- as.data.frame(df3[,2])
colnames(df3)[1] <- "lma_bin"
rownames(df3) <- df_d$name

df4 <- as.data.frame(df_d$ldmc)
colnames(df4)[1] <- "ldmc"
ldmc_quantiles <- unname(ldmc_quantiles)
df4 <- df4 %>% mutate(ldmc_bin = cut(ldmc, breaks=c(min(phy_traits_sp_raw$ldmc) -0.001,ldmc_quantiles[1],ldmc_quantiles[2],
                                                    ldmc_quantiles[3],ldmc_quantiles[4],ldmc_quantiles[5] ,max(phy_traits_sp_raw$ldmc) +0.001)))
df4 <- as.data.frame(df4[,2])
colnames(df4)[1] <- "ldmc_bin"
rownames(df4) <- df_d$name


setdiff(tree3$tip.label, df_d$name)
setdiff(df_d$name, tree3$tip.label)

#plot each trait on phylogeny
circ <- ggtree(tree3, layout = "circular")

p1 <- gheatmap(circ, df1, offset=.8, width= 0.3,
               colnames = F) +
  scale_fill_manual(values = c("(2.86,4.58]"="#E0E9EB",#40004B
                               "(4.58,8.2]"="#ADC5C9",#762A83
                               "(8.2,11.1]"="#7AA3A8",#9970AB
                               "(11.1,16.4]"="#488087",#F1B6DA
                               "(16.4,29]"="#346468",#DE77AE
                               "(29,48.8]"="#21484A"))+
  guides(fill=guide_legend(title=  tau~"(s)"))

p1

p2 <- gheatmap(circ, df2, offset=.8, width= 0.3,
               colnames = F, colnames_angle=0, colnames_offset_y = 150) +
  scale_fill_manual(values = c("(40,42.1]"="#E7D4E8",#fde725
                               "(42.1,44.5]"="#C5A9C8",#7ad151
                               "(44.5,46.1]"="#A47FA9",#22a884
                               "(46.1,47.4]"="#825489",#2a788e
                               "(47.4,49.3]"="#612A6A",#414487
                               "(49.3,51.1]"="#40004B"))+#440154
  guides(fill=guide_legend(title=expression(T[crit]~("°C"))))
p2


p3 <- gheatmap(circ, df3, offset=0, width= 0.3,
               colnames = F, colnames_angle=0, colnames_offset_y = 150) +
  scale_fill_manual(values = c("(0.01,0.019]"="#EEE1D5",#40004B
                               "(0.019,0.031]"="#DEC3AB",#762A83
                               "(0.031,0.048]"="#BD8757",#9970AB
                               "(0.048,0.075]"="#AD692D",#F1B6DA
                               "(0.075,0.144]"="#8C4A21",#DE77AE
                               "(0.144,0.388]"="#6C2C15"))+
  guides(fill=guide_legend(title=bquote("LMA (kg "~m^-2 ~ ")")))

p3

p4 <- gheatmap(circ, df4, offset=0, width= 0.3,
               colnames = F, colnames_angle=0, colnames_offset_y = 150) +
  scale_fill_manual(values = c("(0.033,0.102]"="#EEE1D5",#40004B
                               "(0.102,0.171]"="#DEC3AB",#762A83
                               "(0.171,0.237]"="#BD8757",#9970AB
                               "(0.237,0.293]"="#AD692D",#F1B6DA
                               "(0.293,0.42]"="#8C4A21",#DE77AE
                               "(0.42,0.533]"="#6C2C15"))+
  guides(fill=guide_legend(title=bquote("LDMC (kg "~kg^-1 ~ ")")))

p4

#Save figure
#Note that for tau, the bins within the colour legend are not in the correct order. The order of the bins in the legend is rearranged manually in adobe. 
ridge_plot_supp <- ggarrange(p3, p4, p1, p2,
                             labels = c("a", "b","c", "d"), #hjust = -4, vjust =1, 
                             ncol = 2, nrow = 2)
ridge_plot_supp
ggsave("./figures_revised/Figure3.pdf", plot = ridge_plot_supp, width = 25, height = 20, units = "cm", dpi = 300)

#Statistical tests for differences between angiosperms, gymnosperms and pteridophytes
lma_anova <- aov(lma ~ clade, data = phy_traits_sp)
summary(lma_anova)
TukeyHSD(lma_anova)

ldmc_anova <- aov(ldmc ~ clade, data = phy_traits_sp)
summary(ldmc_anova)

tau_anova <- aov(tau ~ clade, data = phy_traits_sp)
summary(tau_anova)
TukeyHSD(tau_anova)

Tcrit_anova <- aov(Tcrit ~ clade, data = phy_traits_sp)
summary(Tcrit_anova)

Tmax_anova <- aov(Tmax ~ clade, data = phy_traits_sp)
summary(Tmax_anova)
TukeyHSD(Tmax_anova)

#mean and standard error of the mean
mean(phy_traits_sp_raw$lma)
sd(phy_traits_sp_raw$lma)/length(phy_traits_sp_raw$lma)
mean(phy_traits_sp_raw$ldmc)
sd(phy_traits_sp_raw$ldmc)/length(phy_traits_sp_raw$ldmc)
mean(phy_traits_sp_raw$tau)
sd(phy_traits_sp_raw$tau)/length(phy_traits_sp_raw$tau)
mean(phy_traits_sp_raw$Tcrit)
sd(phy_traits_sp_raw$Tcrit)/length(phy_traits_sp_raw$Tcrit)

##### Figure 4: Linear regression #####
#New df to format specifically for phylolm
phy_traits_sp_phylolm<- phy_traits_sp
rownames(phy_traits_sp_phylolm)<-c(phy_traits_sp_phylolm$name)

#PGLS regression for Tcrit and LMA, controlling for variation in mean air temperature on the sampling day
pgls_lma<- phylolm(formula = Tcrit~lma + mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                   boot = 100)

#Null model (no lma)
pgls_nolma<- phylolm(formula = Tcrit~mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                     boot = 100)

#Calculate partial R2 for relationship between Tcrit and LMA 
R2_lik(pgls_lma,pgls_nolma)
summary(pgls_lma)

#Now plot relationship and paste statistical results generated above
corr_lma_Tcrit <- ggplot(phy_traits_sp, aes(x = exp(lma), y = exp(Tcrit))) + #exp() to undo log transformation for plotting absolute values of Tcrit and lma
  
  # Add points for each species
  geom_point( size = 3, alpha = 0.25) +
  #geom_abline(intercept = coefficients["(Intercept)"], slope = coefficients["lma"], color = "red", linetype = "solid") +  # Phylogenetically corrected linear regression line
  # Add phylogenetically corrected linear regression line
  # geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "solid") +
  scale_x_continuous(trans = 'log10')+
  scale_y_continuous(trans = 'log10')+
  
  #geom_point( size = 3, alpha = 0.25) +
  annotate("text", x=0.27, y=50, label= bquote("p = 0.019"), size = 4) + 
    annotate("text", x=0.27, y=51, label= bquote(""~r^2 ~ "= 0.031"), size = 4) + 
  xlab(bquote("LMA (kg "~m^-2 ~ ")")) +
  
  #  xlab("τ (s)") +
  ylab(expression(T[crit]~("°C"))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
corr_lma_Tcrit 

#Now for tau
#PGLS regression for Tcrit and tau, controlling for variation in mean air temperature on the sampling day
pgls_tau<- phylolm(formula = Tcrit~tau + mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                   boot = 100)

#Null model (no tau)
pgls_notau<- phylolm(formula = Tcrit~mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                     boot = 100)

#Calculate partial R2 for relationship between Tcrit and tau 
R2_lik(pgls_tau,pgls_notau)
summary(pgls_tau)


#Now plot relationship and paste statistical results generated above
corr_tau_Tcrit <- ggplot(phy_traits_sp, aes(x = exp(tau), y = exp(Tcrit))) + #exp() to undo log transformation for plotting absolute values of Tcrit and lma
  
  # Add points for each species
  geom_point( size = 3, alpha = 0.25) +
  #geom_abline(intercept = coefficients["(Intercept)"], slope = coefficients["lma"], color = "red", linetype = "solid") +  # Phylogenetically corrected linear regression line
  # Add phylogenetically corrected linear regression line
  # geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "solid") +
  scale_x_continuous(trans = 'log10')+
  scale_y_continuous(trans = 'log10')+
  
  #geom_point( size = 3, alpha = 0.25) +
  annotate("text", x=37, y=50.3, label= bquote("p = 0.351"), size = 4) + 
   annotate("text", x=37, y=51, label= bquote(""~r^2 ~ "= 0.005"), size = 4) + 
  # xlab(bquote("LMA (kg "~m^-2 ~ ")")) +
  
 # xlab("τ (s)") +
  labs(x = tau~"(s)") +
  ylab(expression(T[crit]~("°C"))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
corr_tau_Tcrit

ggarrange(corr_lma_Tcrit,corr_tau_Tcrit,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)

ggsave("./figures_revised/Figure4.pdf", plot = last_plot(), height = 12, width=25 , dpi = 300, units = "cm")


##### Figure 5: PCA ######
numerical_data <- phy_traits_sp[,c(6:12)]
numerical_data <- numerical_data[complete.cases(numerical_data), ]
numerical_data

#rename variables
colnames(numerical_data)[3] <- "LMA"
colnames(numerical_data)[4] <- "LDMC"
colnames(numerical_data)[5] <- "tau"
colnames(numerical_data)[6] <- "L"
colnames(numerical_data)[7] <- "A"

head(numerical_data)

res.pca <- prcomp(numerical_data,  scale = TRUE)
summary(res.pca)
res.pca$rotation

#Define grouping variable for colours
var <- get_pca_var(res.pca)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
#612A6A
#8C4A21

colnames(numerical_data)
#Create PCA plot
pca <- fviz_pca_biplot(res.pca, col.var = colnames(numerical_data), 
                       palette = c("#488087","#488087", "#8C4A21","#8C4A21","#488087", "#612A6A","#612A6A"), #colours labels and arrows according to the category of trait
                       geom = "point",
                       legend.title = "Axis",alpha.ind = 0.5, labelsize = 5, arrowsize = 2)+
  ggtitle("") +
  theme(axis.line = element_line(),
        text = element_text(size = 18),
        panel.grid.major = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
pca

ggsave("./figures_revised/Figure5.pdf", plot = pca,width = 17, height = 15, units = "cm", dpi = 300)



########### Supporting Information ###########
##### Figure S2 ######
covar.1 <- phy_traits_sp[,c(8,9,11)]
head(covar.1)
colnames(covar.1)[1] = "LMA"
colnames(covar.1)[2] = "LDMC"
colnames(covar.1)[3] = "L"

summary(covar.1)

#perform independent effects analysis
IEA <-hier.part(phy_traits_sp$tau, covar.1, family = "gaussian", gof = "Rsqu")
IEA

plotData2_TA <- data.frame(variable = c("LMA", "LDMC","L"), percent = IEA$I.perc$ind.exp.var)
plotData2_TA
plotData2_TA$variable <-factor(plotData2_TA$variable, levels = c("LMA", "LDMC","L"))
plotData2_TA
dev.off()
iea <- ggplot(plotData2_TA)+
  theme(legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ylab("Percentage of independent effects")+
  xlab("")+
  geom_bar(aes(y = percent, x = variable), fill = "#488087", 
           stat = "identity")
iea

ggsave("./figures_revised/Supporting_Info_Figure_S2.pdf", plot = iea, width = 13, height = 10, units = "cm", dpi = 300)

##### Figure S3 ######
#New df to format specifically for phylolm
phy_traits_sp_phylolm<- phy_traits_sp
rownames(phy_traits_sp_phylolm)<-c(phy_traits_sp_phylolm$name)

#PGLS regression for Tcrit and LMA, controlling for variation in mean air temperature on the sampling day
pgls_lma<- phylolm(formula = Tmax~lma + mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                   boot = 100)

#Null model (no lma)
pgls_nolma<- phylolm(formula = Tmax~mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                     boot = 100)

#Calculate partial R2 for relationship between Tcrit and LMA 
R2_lik(pgls_lma,pgls_nolma)
summary(pgls_lma)

corr_lma_Tmax <-ggplot(phy_traits_sp, aes(x = exp(lma), y = exp(Tmax))) + #exp() to undo log transformation for plotting absolute values of Tcrit and lma
  
  # Add points for each species
  geom_point( size = 3, alpha = 0.25) +
  #geom_abline(intercept = coefficients["(Intercept)"], slope = coefficients["lma"], color = "red", linetype = "solid") +  # Phylogenetically corrected linear regression line
  # Add phylogenetically corrected linear regression line
  # geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "solid") +
  scale_x_continuous(trans = 'log10')+
  scale_y_continuous(trans = 'log10')+
  
  #geom_point( size = 3, alpha = 0.25) +
  annotate("text", x=0.02, y=59, label= bquote("p = 0.422"), size = 5) + 
   annotate("text", x=0.02, y=58.2, label= bquote(""~r^2 ~ "= 0.004"), size = 5) + 
  xlab(bquote("LMA (kg "~m^-2 ~ ")")) +
  
  # xlab("τ (s)") +
  ylab(expression(T[max]~("°C"))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
corr_lma_Tmax

### Now for tau

#PGLS regression for Tcrit and tau, controlling for variation in mean air temperature on the sampling day
pgls_tau<- phylolm(formula = Tmax~tau + mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                   boot = 100)

#Null model (no tau)
pgls_notau<- phylolm(formula = Tmax~mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                     boot = 100)

#Calculate partial R2 for relationship between Tcrit and tau 
R2_lik(pgls_tau,pgls_notau)
summary(pgls_tau)

corr_tau_Tmax <- ggplot(phy_traits_sp, aes(x = exp(tau), y = exp(Tmax))) + #exp() to undo log transformation for plotting absolute values of Tcrit and lma
  
  # Add points for each species
  geom_point( size = 3, alpha = 0.25) +
  #geom_abline(intercept = coefficients["(Intercept)"], slope = coefficients["lma"], color = "red", linetype = "solid") +  # Phylogenetically corrected linear regression line
  # Add phylogenetically corrected linear regression line
  # geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "solid") +
  scale_x_continuous(trans = 'log10')+
  scale_y_continuous(trans = 'log10')+
  
  #geom_point( size = 3, alpha = 0.25) +
  annotate("text", x=4.5, y=59, label= bquote("p = 0.053"), size = 5) + 
   annotate("text", x=4.5, y=58.2, label= bquote(""~r^2 ~ "= 0.021"), size = 5) + 
  # xlab(bquote("LMA (kg "~m^-2 ~ ")")) +
  
  labs(x = tau~"(s)") +
  ylab(expression(T[max]~("°C"))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
corr_tau_Tmax

ggarrange(corr_lma_Tmax,corr_tau_Tmax,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)

ggsave("./figures_revised/Supporting_Info_Figure_S3.pdf", plot = last_plot(), height = 12, width=25 , dpi = 300, units = "cm")

##### Figure S4 ######
#PGLS regression for Tcrit and mean temperature on the sampling day
pgls_mean_temp<- phylolm(formula = exp(Tcrit)~mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                         boot = 100)
summary(pgls_mean_temp)
intercept <- as.numeric(coef(pgls_mean_temp)[1])
slope <- as.numeric(coef(pgls_mean_temp)[2])

corr_Tcrit_temp <- ggplot(phy_traits_sp_phylolm, aes(x = mean_temp, y = exp(Tcrit))) + #exp() to undo log transformation for plotting absolute values of Tcrit and lma
  geom_point( size = 3, alpha = 0.25) +
  geom_abline(slope=slope, intercept=intercept, size =2)+
  annotate("text", x=19, y=42, label= bquote("p < 2.2e-16"), size = 5) + 
  annotate("text", x=18.85, y=41.5, label= bquote(""~r^2 ~ "= 0.321"), size = 5) + 
  xlab(expression(Three~day~average~air~temperature~("°C"))) +
  ylab(expression(T[crit]~("°C"))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
corr_Tcrit_temp


#PGLS regression for Tmax and mean temperature on the sampling day
pgls_mean_temp<- phylolm(formula = exp(Tmax)~mean_temp, data = phy_traits_sp_phylolm, phy = tree_up, 
                     boot = 100)
summary(pgls_mean_temp)
intercept <- as.numeric(coef(pgls_mean_temp)[1])
slope <- as.numeric(coef(pgls_mean_temp)[2])
#fitted_df <- data.frame(x = phy_traits_sp_phylolm$mean_temp, y = predict(pgls_mean_temp))

corr_Tmax_temp <- ggplot(phy_traits_sp_phylolm, aes(x = mean_temp, y = exp(Tmax))) + #exp() to undo log transformation for plotting absolute values of Tcrit and lma
  
  geom_point( size = 3, alpha = 0.25) +
  geom_abline(slope=slope, intercept=intercept, size =2)+
  annotate("text", x=19, y=49, label= bquote("p < 2.2e-16 "), size = 5) + 
  annotate("text", x=18.7, y=48.5, label= bquote(""~r^2 ~ "= 0.252"), size = 5) + 
  xlab(expression(Three~day~average~air~temperature~("°C"))) +
  ylab(expression(T[max]~("°C"))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
corr_Tmax_temp

ggarrange(corr_Tcrit_temp,corr_Tmax_temp,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)

ggsave("./figures_revised/Supporting_Info_Figure_S4.pdf", plot = last_plot(), height = 12, width=25 , dpi = 300, units = "cm")


 
##### Figure S5 ######
climate_raw <- read.csv("2022_climate_data.csv")
climate_raw$date <- gsub("-", "", climate_raw$Date.Time)

climate <- climate_raw[,c(8,9,10,12)]

colnames(climate)[1] = "max_temp"
colnames(climate)[2] = "min_temp"
colnames(climate)[3] = "mean_temp"
colnames(climate)[4] = "Date"

results = c()
d= phy_traits_sp
colnames(d)[1] = "Date"

for(j in 1:30) {
  n = j
  x = lag(climate$mean_temp,1)
  len = length(x)
  y = rep(NA, len)
  for (i in 1:len) {
    if (i < n) {
      next
    } else {
      y[i] = mean(x[((i)-(n-1)):(i)], na.rm = T)
    }
  }
  climate$Rolling_Mean_C = y
  clim_cut = climate %>% dplyr::select(Date,Rolling_Mean_C)
  clim_cut$Rolling_Mean_Days = j
  # Match up weather with traits by date
  traits_with_weather = merge(d, clim_cut, by="Date")
  
  # Bind together
  results = bind_rows(results, traits_with_weather)
}

#make 1 column per averaging window
results_wider <- results %>%
  pivot_wider(names_from = Rolling_Mean_Days, values_from = Rolling_Mean_C)


#Calculate R2 for linear regression between Tcrit and mean air temperature by sampling window size
r2_temp_window <- data.frame(matrix(ncol = 3, nrow = 0))
for (i in 1:30){
  window_size = i
  df <- subset(results, results$Rolling_Mean_Days == window_size)
  rownames(df)<-c(df$name)
  lm_Tcrit <- phylolm(formula = Tcrit~Rolling_Mean_C, data = df, phy = tree_up, 
                     boot = 100)
  r2_Tcrit = summary(lm_Tcrit)$adj.r.squared
  lm_Tmax <- phylolm(formula = Tmax~Rolling_Mean_C, data = df, phy = tree_up, 
                      boot = 100)
  r2_Tmax = summary(lm_Tmax)$adj.r.squared
  row <- c(window_size, r2_Tcrit, r2_Tmax)
  r2_temp_window<- rbind(r2_temp_window, row)
}  

x <- c("window_size", "r2_Tcrit", "r2_Tmax")
colnames(r2_temp_window) <- x

r2_bytemp_window_Tcrit <- ggplot(r2_temp_window, aes(x = window_size, y = r2_Tcrit)) +
  geom_line(linewidth = 2, alpha = 0.25)+
  geom_point(size = 3, alpha = 0.75) +  lims(y=c(.13,.32)) +
  xlab("Temperature averaging window (days)") +
  ylab(expression(paste(R^2 ~ (T[crit])))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
r2_bytemp_window_Tcrit


r2_bytemp_window_Tmax <- ggplot(r2_temp_window, aes(x = window_size, y = r2_Tmax)) +
  geom_line(linewidth = 2, alpha = 0.25)+
  geom_point(size = 3, alpha = 0.75) +
  lims(y=c(.13,.32)) +
  xlab("Temperature averaging window (days)") +
  ylab(expression(paste(R^2 ~ (T[max])))) +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'))
r2_bytemp_window_Tmax

ggarrange(r2_bytemp_window_Tcrit,r2_bytemp_window_Tmax,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)

ggsave("./figures_revised/Supporting_Info_Figure_S5.pdf", plot = last_plot(), height = 12, width=25 , dpi = 300, units = "cm")

##### Table S1 ######
#Example code for Tmax, swap trait label to produce moments for other traits
#Taking inverse log as this column was previously log transformed for analysis
min(exp(phy_traits_sp$ldmc), na.rm = T)
min(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "angiosperm"), na.rm = T)
min(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "gymnosperm"), na.rm = T)
min(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "pteridophyte"), na.rm = T)

max(exp(phy_traits_sp$ldmc), na.rm = T)
max(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "angiosperm"), na.rm = T)
max(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "gymnosperm"), na.rm = T)
max(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "pteridophyte"), na.rm = T)

mean(exp(phy_traits_sp$ldmc), na.rm = T)
mean(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "angiosperm"), na.rm = T)
mean(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "gymnosperm"), na.rm = T)
mean(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "pteridophyte"), na.rm = T)

var(exp(phy_traits_sp$ldmc), na.rm = T)
var(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "angiosperm"), na.rm = T)
var(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "gymnosperm"), na.rm = T)
var(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "pteridophyte"), na.rm = T)

skewness(exp(phy_traits_sp$ldmc))
skewness(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "angiosperm"), na.rm = T)
skewness(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "gymnosperm"), na.rm = T)
skewness(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "pteridophyte"), na.rm = T)

kurtosis(exp(phy_traits_sp$ldmc))
kurtosis(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "angiosperm"))
kurtosis(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "gymnosperm"))
kurtosis(subset(exp(phy_traits_sp$ldmc), phy_traits_sp$clade == "pteridophyte"))


###### Table S2 #####
summary(res.pca)

######  Table S3 ####
res.pca$rotation

######  Table S4 ####
#Area
trait <- setNames(phy_traits_sp$leaf_area, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(leaf_area ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)

#Characteristic dimension
trait <- setNames(phy_traits_sp$char_dimension, c(phy_traits_sp$name))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

model<-pgls(char_dimension ~ mean_temp, data=comp.data, lambda = "ML")
summary(model)


######  Table S5 & S6 ####
### Phylo PCA
numerical_data <- phy_traits_sp[,c(6:12)]
numerical_data <- scale(numerical_data) 
rownames(numerical_data) <- phy_traits_sp$name
phy_pca <- phyl.pca(tree3, numerical_data, method="BM", mode="cov")

## S3 method for class 'phyl.pca'
biplot(phy_pca )
summary(phy_pca)
phy_pca$L
