
        ##################
        #manuscript plots#
        ##################
        
#we do not want to include prs noise info

library(ggplot2)
library(patchwork)
library(ggpubr)

##############################
# regression model III plots #
##############################
        
reg <- read.table("sim_noIGE_trioprs.txt", header=T)
colnames(reg) <- c("Estimate", "SE", "T-value", "P-value", "AdjR2", "prs_noise", "replica", "h2", "Tagged")
reg <- subset(reg, reg$prs_noise==0)

df <- reg
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)

#maternal

mat <- seq(3, 540, by=4)
df1 <- df[mat,]

plot1 <- ggplot(data=df1, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Model III; Maternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=20), 
        strip.text = element_text(size=20), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm3_mIGE.png", plot=plot1, width=10, height=5)
ggsave(file="regm3_mIGE.svg", plot=plot1, width=10, height=5)

#paternal

pat <- seq(4, 540, by=4)
df2 <- df[pat,]

plot2 <- ggplot(data=df2, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Model III; Paternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="regm3_pIGE.png", plot=plot2, width=10, height=5)
ggsave(file="regm3_pIGE.svg", plot=plot2, width=10, height=5)

#############################
# regression model IV plots #
#############################

reg <- read.table("sim_noIGE_trioprs_parentpheno.txt", header=T)
colnames(reg) <- c("Predictor", "Estimate", "SE", "T-value", "P-value", "AdjR2", "prs_noise", "replica", "h2", "Tagged")
reg <- subset(reg, reg$prs_noise==0)

df <- reg
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)

#maternal (1)

df3 <- df[df$Predictor=="maternal_prs",]

plot3 <- ggplot(data=df3, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Maternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="regm4_mIGE.png", plot=plot3, width=10, height=5)
ggsave(file="regm4_mIGE.svg", plot=plot3, width=10, height=5)

#paternal (1)

df4 <- df[df$Predictor=="paternal_prs",]

plot4 <- ggplot(data=df4, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Paternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="regm4_pIGE.png", plot=plot4, width=10, height=5)
ggsave(file="regm4_pIGE.svg", plot=plot4, width=10, height=5)

#maternal (2)

df5 <- df[df$Predictor=="maternal_phenotype",]

plot5 <- ggplot(data=df5, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Maternal Phenotype") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="regm4_mIGE_2.png", plot=plot5, width=10, height=5)
ggsave(file="regm4_mIGE_2.svg", plot=plot5, width=10, height=5)

#paternal (2)

df6 <- df[df$Predictor=="paternal_phenotype",]

plot6 <- ggplot(data=df6, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Paternal Phenotype") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="regm4_pIGE_2.png", plot=plot6, width=10, height=5)
ggsave(file="regm4_pIGE_2.svg", plot=plot6, width=10, height=5)

#combine maternal plots

mplots <- plot1 + plot3 + plot5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

mplots[[1]] <- mplots[[1]] + theme(axis.title.x = element_blank())

mplots[[2]] <- mplots[[2]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank())

mplots[[3]] <- mplots[[3]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank())

#annotating figure to include main title

mplots <- mplots + plot_annotation(title = "Maternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.45)))


ggsave(file="maternal_rGE_reg_combined_manuscript.png", plot=mplots, width=15, height=6)
ggsave(file="maternal_rGE_reg_combined_manuscript.svg", plot=mplots, width=15, height=6)

#combine paternal plots

pplots <- plot2 + plot4 + plot6 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

pplots[[1]] <- pplots[[1]] + theme(axis.title.x = element_blank())

pplots[[2]] <- pplots[[2]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank())

pplots[[3]] <- pplots[[3]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank())

#annotating figure to include main title

pplots <- pplots + plot_annotation(title = "Paternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.45)))


ggsave(file="paternal_rGE_reg_combined_manuscript.png", plot=pplots, width=15, height=6)
ggsave(file="paternal_rGE_reg_combined_manuscript.svg", plot=pplots, width=15, height=6)

######################
#pathway model plots #
######################

######################
# simple model plots #
######################

m1 <- read.table("sim_pa_noIGE_m1_relabelled.txt", header=T)
colnames(m1) <- c("DV", "IV", "label", "Estimate", "SE", "P-value", "lowerCI", "upperCI", "AIC", "replica", "h2", "Tagged", "prs_noise")
m1 <- subset(m1, m1$prs_noise==0)

df <- m1
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)

#paternal genetic nurturing

df1 <- df[df$label=="b",]

plot1 <- ggplot(data=df1, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Simple Model; Path b") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="sim_simple_pIGE_manuscript.png", plot=plot1, width=10, height=5)
ggsave(file="sim_simple_pIGE_manuscript.svg", plot=plot1, width=10, height=5)

#maternal genetic nurturing

df2 <- df[df$label=="c",]

plot2 <- ggplot(data=df2, aes(x=h2, y=Estimate, fill=Tagged)) + 
  geom_boxplot(outlier.size=0, fatten=0.5) + 
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Beta Estimate") +
  ggtitle("Simple Model; Path c") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="sim_simple_mIGE_manuscript.png", plot=plot2, width=10, height=5)
ggsave(file="sim_simple_mIGE_manuscript.svg", plot=plot2, width=10, height=5)

########################
# extended model plots #
########################

m2 <- read.table("sim_pa_noIGE_m2_relabelled.txt", header=T)
colnames(m2) <- c("DV", "IV", "label", "Estimate", "SE", "P-value", "lowerCI", "upperCI", "AIC", "replica", "h2", "Tagged", "prs_noise")
m2 <- subset(m2, m2$prs_noise==0)

df <- m2
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)

#paternal genetic nurturing (1)

df3 <- df[df$label=="b",]

plot3 <- ggplot(data=df3, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Offspring Pheno ~ Paternal PGS Beta Estimate") +
  ggtitle("Extended Model; Path b") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="sim_extended_pIGE_manuscript.png", plot=plot3, width=10, height=5)
ggsave(file="sim_extended_pIGE_manuscript.svg", plot=plot3, width=10, height=5)

#maternal genetic nurturing (1)

df4 <- df[df$label=="c",]

plot4 <- ggplot(data=df4, aes(x=h2, y=Estimate, fill=Tagged)) + 
  geom_boxplot(outlier.size=0, fatten=0.5) + 
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Offspring Pheno ~ Maternal PGS Beta Estimate") +
  ggtitle("Extended Model; Path c") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="sim_extended_mIGE_manuscript.png", plot=plot4, width=10, height=5)
ggsave(file="sim_extended_pIGE_manuscript.svg", plot=plot4, width=10, height=5)

#paternal genetic nurturing (2)

df5 <- df[df$label=="f",]

plot5 <- ggplot(data=df5, aes(x=h2, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Offspring Pheno ~ Paternal Pheno Estimate") +
  ggtitle("Extended Model; Path f") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="sim_extended_pIGE2_manuscript.png", plot=plot5, width=10, height=5)
ggsave(file="sim_extended_pIGE2_manuscript.svg", plot=plot5, width=10, height=5)

#maternal genetic nurturing (2)

df6 <- df[df$label=="g",]

plot6 <- ggplot(data=df6, aes(x=h2, y=Estimate, fill=Tagged)) + 
  geom_boxplot(outlier.size=0, fatten=0.5) + 
  scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Trait Heritability") + 
  ylab("Offspring Pheno ~ Maternal Pheno Beta Estimate") +
  ggtitle("Extended Model; Path g") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), strip.text = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=20))

ggsave(file="sim_extended_mIGE2_manuscript.png", plot=plot6, width=10, height=5)
ggsave(file="sim_extended_mIGE2_manuscript.svg", plot=plot6, width=10, height=5)

#maternal plots#
#combine maternal genetic nurturing path plots and creating a common legend

mplots <- plot2 + plot4 + plot6 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

mplots[[1]] <- mplots[[1]] + theme(axis.title.x = element_blank())

mplots[[2]] <- mplots[[2]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank())

mplots[[3]] <- mplots[[3]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank())

#annotating figure to include main title

mplots <- mplots + plot_annotation(title = "Maternal Genetic Nurturing Effects",
                         theme = theme(plot.title = element_text(size = 20, hjust = 0.50)))


ggsave(file="maternal_rGE_pa_combined_manuscript.png", plot=mplots, width=15, height=6)
ggsave(file="maternal_rGE_pa_combined_manuscript.svg", plot=mplots, width=15, height=6)

#paternal plots#
#combine paternal genetic nurturing path plots and creating a common legend

pplots <- plot1 + plot3 + plot5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

pplots[[1]] <- pplots[[1]] + theme(axis.title.x = element_blank())

pplots[[2]] <- pplots[[2]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank())

pplots[[3]] <- pplots[[3]] + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.title.x = element_blank())

#annotating figure to include main title

pplots <- pplots + plot_annotation(title = "Paternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.50)))


ggsave(file="paternal_rGE_pa_combined_manuscript.png", plot=pplots, width=15, height=6)
ggsave(file="paternal_rGE_pa_combined_manuscript.svg", plot=pplots, width=15, height=6)

