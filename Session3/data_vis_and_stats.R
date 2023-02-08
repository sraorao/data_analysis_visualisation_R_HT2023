###R for biologists
##Irina & Rao, 22/02/2023
#Import and prepare the data ####
library(tidyverse)
#Import the all sheets from Excel file provided (use readxl)
library(readxl)
ifny <- read_excel("Session2/data/MSD_data.xlsx", sheet = 'IFNy')
il1 <- read_excel("Session2/data/MSD_data.xlsx", sheet = 'IL-1')
il2 <- read_excel("Session2/data/MSD_data.xlsx", sheet = 'IL-2')
il10 <- read_excel("Session2/data/MSD_data.xlsx", sheet = 'IL-10')
library(dplyr)
ifny %>% mutate(Subtraction=Calc.Conc.Mean_S1S2-Calc.Conc.Mean_Unst) -> ifny
il1 %>% mutate(Subtraction=Calc.Conc.Mean_S1S2-Calc.Conc.Mean_Unst) -> il1
il2 %>% mutate(Subtraction=Calc.Conc.Mean_S1S2-Calc.Conc.Mean_Unst) -> il2
il10 %>% mutate(Subtraction=Calc.Conc.Mean_S1S2-Calc.Conc.Mean_Unst) -> il10

ifny<-ifny[,c(1,5,6)]
il1<-il1[,c(1,5,6)]
il2<-il2[,c(1,5,6)]
il10<-il10[,c(1,5,6)]

cyto <- Reduce(function(...) merge(..., all = TRUE, by = c("Participant_ID","Randomised_to")),
       list(ifny, il1, il2, il10))
colnames(cyto)<-c("ID","Randomisation","IFNy","IL1","IL2","IL10")
#Remove negative values
cyto %>% filter(across(c(IFNy, IL1, IL2, IL10),~ .x>0|is.na(.x)==TRUE)) -> cyto
cyto %>% filter(across(c(-ID,-Randomisation),~ .x>0|is.na(.x)==TRUE)) -> cyto

#Convert groups into factors
cyto$Randomisation <- as.factor(cyto$Randomisation)
levels(cyto$Randomisation)<-c("ChAdOx1", "Control")

cyto_melted<-gather(cyto,key="cytokine",value = "value",-ID,-Randomisation)
cyto_melted<-na.omit(cyto_melted)

#Continue plotting from last session ####
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)
#Dot+boxplots
ggplot(cyto_melted, aes(fill=Randomisation, y=log(value), x=cytokine)) +
  geom_boxplot(width = 0.8, outlier.shape = NA) + 
  geom_dotplot(binaxis = "y", stackdir = "center",
               position = position_dodge(0.8),dotsize = 0.8, binwidth = 0.2) + 
  theme_pubr() #publication ready theme

#Jitter plots
ggplot(cyto_melted, aes(fill=Randomisation, y=log(value), x=cytokine))+
geom_jitter(aes(shape=Randomisation, color=Randomisation),
  position = position_jitter(0.2),
  size = 1.2)+
  theme_pubr()

#Change positions 
ggplot(cyto_melted, aes(fill=Randomisation, y=log(value), x=cytokine))+
  geom_jitter(aes(shape=Randomisation, color=Randomisation),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.2)+
  theme_pubr() 

#Adding stat
ggplot(cyto_melted, aes(fill=Randomisation, y=log(value), x=cytokine)) +
  geom_jitter(aes(shape=Randomisation, color=Randomisation),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.2) + 
  theme_pubr() + #publication-ready theme
  stat_summary(
  fun.data="mean_sdl",  fun.args = list(mult=1), 
  geom = "pointrange",  size = 0.4,
  position = position_dodge(0.8), show.legend = FALSE)

#Finalising graph format
ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2) + 
  theme_pubr() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.5), show.legend = FALSE)


#Performing statistics and adding significance levels####
#Summary
compare_means(data = cyto_melted, formula = value ~ Randomisation, group.by = "cytokine")
#Adding to the ggplot
#Boxplot
ggplot(cyto_melted, aes(fill=Randomisation, y=log(value), x=cytokine)) + 
  geom_boxplot() +
  theme_pubr() +
  stat_compare_means(aes(group = Randomisation), method = "wilcox", hide.ns = F, paired = F,
                     label = "p.signif",bracket.size = 0.3)
#alternatively - label = "p.format"

#Barplot
ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  geom_bar(position="dodge", stat = "summary", fun = "mean") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, 
               position = position_dodge(width = 0.90), width = 0.3) +
  stat_compare_means(aes(group = Randomisation), label = "p.signif", label.y = 75)+
  theme_bw()


library(rstatix)
#Use p-adjusted values instead
stat.test <- cyto_melted %>%
  group_by(cytokine) %>%
  wilcox_test(value ~ Randomisation) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
View(stat.test)
stat.test <- stat.test %>% add_xy_position(x = "cytokine")

ggplot(cyto_melted, aes(fill=Randomisation, y=value, x=cytokine)) + 
  geom_bar(position="dodge", stat = "summary", fun = "mean") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, 
               position = position_dodge(width = 0.90), width = 0.3) +
  theme_bw() + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0,hide.ns = TRUE, inherit.aes = FALSE, y.position = 80)


#Finalising graph format + stats
ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2) + 
  theme_pubr() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.5), show.legend = FALSE) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0,hide.ns = TRUE, inherit.aes = FALSE, y.position = 3)

#Create interactive plot
# install.packages("plotly")
library(plotly)
p<-ggplot(cyto_melted, aes(group=Randomisation, y=value, x=cytokine)) +
  scale_y_log10() +
  geom_jitter(aes(color=Randomisation),
              position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
              size = 1.2) + 
  theme_pubr() +
  scale_color_manual(values=c("red","blue")) +
  xlab("") +
  ylab("Cytokine concentration") +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.5), show.legend = FALSE) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0,hide.ns = TRUE, inherit.aes = FALSE, y.position = 3)

ggplotly(p)
