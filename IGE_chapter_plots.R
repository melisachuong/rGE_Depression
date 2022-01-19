


            ##################
            # chapter plots  #
            ##################

            
library(ggplot2)
library(patchwork)
library(ggpubr)

#############################            
# regression model III plots#
#############################
            
reg <- read.table("sim_noIGE_trioprs.txt", header=T)
colnames(reg) <- c("Estimate", "SE", "T-value", "P-value", "AdjR2", "prs_noise", "replica", "h2", "Tagged")

df <- reg
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)
df$prs_noise <- as.factor(df$prs_noise)

#maternal

mat <- seq(3, 5940, by=4)
df1 <- df[mat,]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot1 <- ggplot(data=df1, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN,limits=c(-0.15,0.15), breaks=seq(-0.15,0.15, 0.075), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Model III; Maternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm3_mIGE_chapter.png", plot=plot1, width=15, height=5)
ggsave(file="regm3_mIGE_chapter.svg", plot=plot1, width=15, height=5)
#regression model 3, maternal IGE

#paternal

pat <- seq(4, 5940, by=4)
df2 <- df[pat,]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot2 <- ggplot(data=df2, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN,limits=c(-0.15,0.15), breaks=seq(-0.15,0.15, 0.05), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Model III; Paternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm3_pIGE_chapter.png", plot=plot2, width=15, height=5)
ggsave(file="regm3_pIGE_chapter.svg", plot=plot2, width=15, height=5)

#############################
# regression model IV plots #
#############################

reg <- read.table("sim_noIGE_trioprs_parentpheno.txt", header=T)
colnames(reg) <- c("Predictor", "Estimate", "SE", "T-value", "P-value", "AdjR2", "prs_noise", "replica", "h2", "Tagged")

df <- reg
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)
df$prs_noise <- as.factor(df$prs_noise)


#maternal (1)

df3 <- df[df$Predictor=="maternal_prs",]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot3 <- ggplot(data=df3, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN,limits=c(-0.5,0.1), breaks=seq(-0.5,0.1, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Maternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm4_mIGE_chapter.png", plot=plot3, width=15, height=5)
ggsave(file="regm4_mIGE_chapter.svg", plot=plot3, width=10, height=5)

#paternal (1)

df4 <- df[df$Predictor=="paternal_prs",]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot4 <- ggplot(data=df4, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN,limits=c(-0.5, 0.1), breaks=seq(-0.5,0.1, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Paternal Polygenic Score") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm4_pIGE_chapter.png", plot=plot4, width=15, height=5)
ggsave(file="regm4_pIGE_chapter.svg", plot=plot4, width=10, height=5)

#maternal (2)

df5 <- df[df$Predictor=="maternal_phenotype",]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot5 <- ggplot(data=df5, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN,limits=c(-0.1,0.5), breaks=seq(-0.1,0.5, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Maternal Phenotype") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm4_mIGE_2_chapter.png", plot=plot5, width=15, height=5)
ggsave(file="regm4_mIGE_2_chapter.svg", plot=plot5, width=10, height=5)

#paternal (2)

df6 <- df[df$Predictor=="paternal_phenotype",]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot6 <- ggplot(data=df6, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN,limits=c(-0.1,0.5), breaks=seq(-0.1,0.5, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Model IV; Paternal Phenotype") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="regm4_pIGE_2_chapter.png", plot=plot5, width=15, height=5)
ggsave(file="regm4_pIGE_2_chapter.svg", plot=plot5, width=10, height=5)

#maternal plots#
#combine maternal genetic nurturing path plots and creating a common legend

mplots <- plot1 / plot3 
mplots <- mplots/ plot5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

mplots[[1]] <- mplots[[1]] + theme(axis.title.y = element_blank(),
                                   axis.title.x = element_blank())

mplots[[2]] <- mplots[[2]] + theme(axis.title.x = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())

mplots[[3]] <- mplots[[3]] + theme(axis.title.y = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())
#annotating figure to include main title

mplots <- mplots + plot_annotation(title = "Regression Models; Maternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.50)))

ggsave(file="maternal_rGE_reg_combined_chapter.png", plot=mplots, width=20, height=20)
ggsave(file="maternal_rGE_reg_combined_chapter.svg", plot=mplots, width=20, height=20)

#paternal plots#
#combine paternal genetic nurturing path plots and creating a common legend

pplots <- plot2 / plot4 
pplots <- pplots/ plot6 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

pplots[[1]] <- pplots[[1]] + theme(axis.title.y = element_blank(),
                                   axis.title.x = element_blank())

pplots[[2]] <- pplots[[2]] + theme(axis.title.x = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())

pplots[[3]] <- pplots[[3]] + theme(axis.title.y = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())
#annotating figure to include main title

pplots <- pplots + plot_annotation(title = "Regression Models; Paternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.50)))

ggsave(file="paternal_rGE_reg_combined_chapter.png", plot=pplots, width=20, height=20)
ggsave(file="paternal_rGE_reg_combined_chapter.svg", plot=pplots, width=20, height=20)



        #####################
        #pathway model plots#
        #####################

######################
# simple model plots #
######################
m1 <- read.table("sim_pa_noIGE_m1_relabelled.txt", header=T)
colnames(m1) <- c("DV", "IV", "label", "Estimate", "SE", "P-value", "lowerCI", "upperCI", "AIC", "replica", "h2", "Tagged", "prs_noise")

df <- m1
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)
df$prs_noise <- as.factor(df$prs_noise)

#paternal genetic nurturing

df1 <- df[df$label=="b",]

require(scales)
scaleFUN <- function(x) sprintf("%.2f", x)

plot1 <- ggplot(data=df1, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(labels = scaleFUN, limits=c(-0.15,0.15), breaks=seq(-0.15,0.15, 0.075), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Simple Model; Path b") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) +
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=20), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="sim_simple_pIGE_chapter.png", plot=plot1, width=10, height=5)
ggsave(file="sim_simple_pIGE_chapter.svg", plot=plot1, width=10, height=5)

#maternal genetic nurturing

df2 <- df[df$label=="c",]

plot2 <- ggplot(data=df2, aes(x=prs_noise, y=Estimate, fill=Tagged)) + 
  geom_boxplot(outlier.size=0, fatten=0.5) + 
  scale_y_continuous(labels = scaleFUN, limits=c(-0.15,0.15), breaks=seq(-0.15,0.15, 0.075), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Simple Model; Path c") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) +
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="sim_simple_mIGE_chapter.png", plot=plot1, width=10, height=5)
ggsave(file="sim_simple_mIGE_chapter.svg", plot=plot1, width=10, height=5)

########################
# extended model plots #
########################

m2 <- read.table("sim_pa_noIGE_m2_relabelled.txt", header=T)
colnames(m2) <- c("DV", "IV", "label", "Estimate", "SE", "P-value", "lowerCI", "upperCI", "AIC", "replica", "h2", "Tagged", "prs_noise")

df <- m2
df$Tagged <- as.factor(df$Tagged)
df$h2 <- as.factor(df$h2)
df$prs_noise <- as.factor(df$prs_noise)

#paternal genetic nurturing (1)

df3 <- df[df$label=="b",]

plot3 <- ggplot(data=df3, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.5,0.1), breaks=seq(-0.5,0.1, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Extended Model; Path b") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) +
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=20), 
        strip.text = element_text(size=20), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))


ggsave(file="sim_extended_pIGE_chapter.png", plot=plot3, width=10, height=5)
ggsave(file="sim_extended_pIGE_chapter.svg", plot=plot3, width=10, height=5)

#maternal genetic nurturing (1)

df4 <- df[df$label=="c",]

plot4 <- ggplot(data=df4, aes(x=prs_noise, y=Estimate, fill=Tagged)) + 
  geom_boxplot(outlier.size=0, fatten=0.5) + 
  scale_y_continuous(limits=c(-0.5,0.1), breaks=seq(-0.5,0.5, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Extended Model; Path c") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) +
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="sim_extended_mIGE_chapter.png", plot=plot3, width=10, height=5)
ggsave(file="sim_extended_mIGE_chapter.svg", plot=plot3, width=10, height=5)

#paternal genetic nurturing (2)

df5 <- df[df$label=="f",]

plot5 <- ggplot(data=df5, aes(x=prs_noise, y=Estimate, fill=Tagged)) +
  geom_boxplot(outlier.size=0, fatten=0.5) +
  scale_y_continuous(limits=c(-0.1,0.5), breaks=seq(-0.1,0.5, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Extended Model; Path f") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, size=20), 
        strip.text = element_text(size=20), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="sim_extended_pIGE2_chapter.png", plot=plot5, width=10, height=5)
ggsave(file="sim_extended_pIGE2_chapter.svg", plot=plot5, width=10, height=5)

#maternal genetic nurturing (2)

df6 <- df[df$label=="g",]

plot6 <- ggplot(data=df6, aes(x=prs_noise, y=Estimate, fill=Tagged)) + 
  geom_boxplot(outlier.size=0, fatten=0.5) + 
  scale_y_continuous(limits=c(-0.1,0.5), breaks=seq(-0.1,0.5, 0.1), expand = c(0, 0)) +
  xlab("Polygenic Score Noise") + 
  ylab("Beta Estimate") +
  ggtitle("Extended Model; Path g") +
  scale_fill_brewer("Tagged Genetic Variance", palette="Accent") + 
  facet_wrap(~h2) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20), 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, size=15), 
        strip.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

ggsave(file="sim_extended_mIGE2.png", plot=plot6, width=10, height=5)
ggsave(file="sim_extended_mIGE2.svg", plot=plot6, width=10, height=5)

#maternal plots#
#combine maternal genetic nurturing path plots and creating a common legend

mplots <- plot2 / plot4

mplots <- mplots / plot6 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

mplots[[1]] <- mplots[[1]] + theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank())

mplots[[2]] <- mplots[[2]] + theme(axis.title.x = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())

mplots[[3]] <- mplots[[3]] + theme(axis.title.y = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())

#annotating figure to include main title

mplots <- mplots + plot_annotation(title = "Maternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.50)))

ggsave(file="maternal_rGE_pa_combined_chapter.png", plot=mplots, width=15, height=15)
ggsave(file="maternal_rGE_pa_combined_chapter.svg", plot=mplots, width=15, height=15)


#paternal plots#
#combine paternal genetic nurturing path plots and creating a common legend

pplots <- plot1 / plot3

pplots <- pplots / plot5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

#create a common y-axis and removing axis-titles

pplots[[1]] <- pplots[[1]] + theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank())

pplots[[2]] <- pplots[[2]] + theme(axis.title.x = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())

pplots[[3]] <- pplots[[3]] + theme(axis.title.y = element_blank(),
                                   strip.background = element_blank(),
                                   strip.text.x = element_blank())

#annotating figure to include main title

pplots <- pplots + plot_annotation(title = "Paternal Genetic Nurturing Effects",
                                   theme = theme(plot.title = element_text(size = 20, hjust = 0.50)))

ggsave(file="paternal_rGE_pa_combined_chapter.png", plot=pplots, width=15, height=15)
ggsave(file="paternal_rGE_pa_combined_chapter.svg", plot=pplots, width=15, height=15)

