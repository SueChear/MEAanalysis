library(magrittr)
library(stringr)
library(dplyr)
library(ggplot2)


setwd("P:/PhD project/53. final MEA/extracted csv 2")

#to read the list of doseresponse csv, then extract the DIV and batch string characters
#in the file title, then cbind new column "DIV" into the main csv file


my_read_csv <- function(x) {
  out <- read.csv(x,)
  div <- gsub(".*DIV(.+)_B.*","\\1",x)
  batch<-gsub(".*_B(.+).csv","\\1",x)
  cbind(DIV=div, Batch=batch, out)
}

#list the dose response files
db_list <- list.files("P:/PhD project/53. final MEA/extracted csv 2", recursive = T, pattern="doseresponse", full.names = F)
tbl <- lapply(db_list, my_read_csv) %>% dplyr::bind_rows()
data<-tbl
View(data)

#replace all NA values with 0
data[is.na(data)] <- 0

#remove certain irregular characters in colnames
colnames(data)<-gsub("..µs.","",colnames(data))

colnames(data)<-gsub("?..","",colnames(data))
colnames(data)

#select only active channels
data2<-subset(data, Active.Channel=="True")
View(data2)


#change burst count to burst count/min, and burst duration to sec
data3<-data2%>%dplyr::mutate(Burst.Rate=Burst.Count/5,
                             Mean.Burst.Duration=Mean.Burst.Duration/1000000,
                             Mean.Interburst.Interval=Mean.Interburst.Interval/1000000,
                             Mean.Network.Burst.Duration=Mean.Network.Burst.Duration/1000000,
                             Mean.Network.Interburst.Interval=Mean.Network.Interburst.Interval/1000000,
                             Network.Burst.Rate=Network.Burst.Count/5)
View(data3)                                 


data4<-mutate(data3,
              Well_Batch = paste(Well.Label, Batch,sep = '_'))

View(data4)
#remove unwanted columns 
data5<-subset(data4, select=-c(Channel.Label,Batch,Well.ID,Well.Label,Compound.Name,Experiment,
                               Dose..pM.,Dose.Label,Spike.Count,Active.Channel,Network.Burst.Count,Burst.Count,
                               Mean.Burst.Spike.Count,Mean.Network.Burst.Spike.Count))
View(data5)

#summarize into averages of various parameters
colnames(data5)
data6<-data5%>% dplyr::group_by(DIV,Well_Batch,Compound.ID) %>% 
  summarize(Electrode.n = n(),Mean.firing.rate = mean(Spike.Rate..Hz.),Mean.burst.duration=mean(Mean.Burst.Duration),
            Mean.Burst.Spike.Rate=mean(Mean.Burst.Spike.Rate..Hz.),Percentage.Spikes.in.Burst=mean(Percentage.Spikes.in.Burst),
            Mean.Interburst.Interval=mean(Mean.Interburst.Interval),Mean.Network.Burst.Duration=mean(Mean.Network.Burst.Duration),
            Mean.Network.Burst.Spike.Rate=mean(Mean.Network.Burst.Spike.Rate..Hz.),
            Percentage.Spikes.in.Network.Burst=mean(Percentage.Spikes.in.Network.Burst),
            Mean.Network.Interburst.Interval=mean(Mean.Network.Interburst.Interval),
            Mean.burst.rate=mean(Burst.Rate), Mean.Network.Burst.Rate=mean(Network.Burst.Rate))%>%
  arrange(DIV,Compound.ID)
View(data6)


#arrange DIV into the order below

x <- c("4", "5", "","6","7","8","9","10","11","12","13","14","15","16","17",
       "18","19","20","21","22","23","24","25","26","27","28","29","30",
       "31","32","33","34","35","36","37","38","39","40","41","42")

data7<-data6 %>%dplyr::mutate(DIV =  factor(DIV, levels = x)) %>%
  arrange(DIV)    
View(data7)

data7<-data7%>%filter(Compound.ID!="TOB220")
View(data7)


#function for timestamp csvs, read the list, extract the title DIV
my_read_csv <- function(x) {
  out <- read.csv(x,)
  div <- gsub(".*DIV(.+)_B.*","\\1",x)
  batch<-gsub(".*_B(.+).csv","\\1",x)
  cbind(DIV=div, Batch=batch, out)
}

db_list2 <- list.files("P:/PhD project/53. final MEA/extracted csv 2", recursive = T, pattern="timestamp", full.names = F)
tbl2 <- lapply(db_list2, my_read_csv) %>% dplyr::bind_rows()
dataT<-tbl2

View(dataT)
dataT2<-dataT%>%dplyr::rename(Timestamp=Timestamp..µs.,
                              Channel.ID=?..Channel.ID)
View(dataT2)

dataT3<-mutate(dataT2,
               Well_Batch = paste(Well.Label, Batch,sep = '_'))
View(dataT3)

#select columns
dataT4<-dataT3%>%select(DIV,Well_Batch,Channel.ID,Compound.ID,Timestamp)
View(dataT4)

#deduct timestamp n from timestamp n+1
dataT5<-dataT4%>%group_by(DIV,Well_Batch,Channel.ID) %>%
  mutate(diff = Timestamp - lag(Timestamp, default = first(Timestamp)))
View(dataT5)

#remove initial timestamp interval 0
dataT6<-dataT5[dataT5$diff != 0, ]
View(dataT6)

#merge data3 per electrode data with ISI
#average the intervals
dataT6.2<-dataT6%>% group_by(DIV,Well_Batch,Channel.ID,Compound.ID) %>% summarize(Mean.ISI = mean(diff))%>%
  arrange(DIV,Compound.ID)
View(dataT6.2)

#re-arrange DIV to the order specified by x
x <- c("4", "5", "","6","7","8","9","10","11","12","13","14","15","16","17",
       "18","19","20","21","22","23","24","25","26","27","28","29","30",
       "31","32","33","34","35","36","37","38","39","40","41","42")
dataT7.2<-dataT6.2 %>%dplyr::mutate(DIV =  factor(DIV, levels = x)) %>%
  arrange(DIV,Compound.ID) 
View(dataT7.2)

dataT7.2<-dataT7.2%>%filter(Compound.ID!="TOB220")



#list out active wells 
list_wells<-data.frame(data5$DIV,data5$Well_Batch,data5$Compound.ID,data5$Channel.ID)
colnames(list_wells) <- c("DIV","Well_Batch","Compound.ID","Channel.ID")
View(list_wells)

#merge active wells to ISI
joined_data2 <- dataT7.2 %>% inner_join(list_wells, by=c("DIV","Well_Batch","Channel.ID"))
View(joined_data2)     

#remove not needed columns
joined_data2 <- joined_data2[ -c(6) ]

#rename columns
joined_data2<-joined_data2%>%dplyr::rename(Compound.ID=Compound.ID.x)
View(joined_data2)

#average mean ISI over channel number

joined_data3<-joined_data2%>% dplyr::group_by(DIV,Well_Batch,Compound.ID) %>% 
  summarize(Mean.ISI = mean(Mean.ISI))%>%
  arrange(DIV,Compound.ID)
View(joined_data3)


#join ISI to main dataset


data_join2<-data7 %>% inner_join(joined_data3, by=c("DIV","Well_Batch","Compound.ID"))
View(data_join2)
data_join2.2<-dplyr::mutate(data_join2, MISI=Mean.ISI/1000000)
View(data_join2.2)

#add column for percentage of active electrodes
data_join3<-data_join2.2%>%dplyr::mutate(Percent.Act.Elec=100*Electrode.n/12)
View(data_join3)

#remove not needed columns
data_join3 <- data_join3[ -c(16) ]

#rename column
colnames(data_join3)
data_final<-data_join3%>%dplyr::rename(Cell.line=Compound.ID,
                                       'Number of active electrode'=Electrode.n,
                                       'Mean.Firing.Rate(Hz)'=Mean.firing.rate,
                                       'Mean.Burst.Duration(s)'=Mean.burst.duration,
                                       'Mean.Burst.Spike.Rate(Spikes/s)'=Mean.Burst.Spike.Rate,
                                       'Mean.Interburst.Interval(s)'=Mean.Interburst.Interval,
                                       'Mean.Network.Burst.Duration(s)'=Mean.Network.Burst.Duration,
                                       'Mean.Network.Burst.Spike.Rate(Spikes/s)'=Mean.Network.Burst.Spike.Rate,
                                       'Mean.Network.Interburst.Interval(s)'=Mean.Network.Interburst.Interval,
                                       'Mean.burst.rate(Burst/min)'=Mean.burst.rate,
                                       'Mean.Network.Burst.Rate(Network burst/min)'=Mean.Network.Burst.Rate,
                                       'Percentage of active electrode'=Percent.Act.Elec,
                                       'Mean.ISI(s)'=MISI)
View(data_final)
#replace dots with space before exporting csv
data_exp<-data_final
names(data_exp) <- gsub(x = names(data_exp),
                        pattern = "\\.",
                        replacement = " ")
#final dataframe
View(data_exp)

#write csv
write.csv(data_exp,"data_final_2.csv", row.names = FALSE)

####start here if csv available####
setwd("P:/PhD project/71. organized thesis data/MEA")
sink(file="MEA_2_R.txt")
data_final_2<-read.csv("data_final_2.csv", header=T, sep=",")


#analysis####
## fit Generalized Additive Mixed-effect Models (GAMM)

# location-scale model assuming gamma-distributed residuals
library(mgcv)
library(plyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(magrittr)

data_final_2$Cell.line <- factor(data_final_2$Cell.line)
data_final_2$Well_Batch <- factor(data_final_2$Well_Batch)




cols <- names(data_final_2)[5:ncol(data_final_2)]
View(cols)
result1 <- lapply(cols, function(z) mgcv::gam(list(as.formula(paste0(z,"+1~s(DIV, by = Cell.line, bs = 'cs') + s(Well_Batch, bs = 're')")),
                                                   ~s(DIV, by = Cell.line, bs = 'cs')),
                                              data = data_final,
                                              family = gammals()))

# assume gamma-distributed residuals
result <- lapply(cols, function(z) mgcv::gam(as.formula(paste0(z,"+1~s(DIV, by = Cell.line, bs = 'cs') + s(Well_Batch, bs = 're')")),
                                             data = data_final_2,
                                             family = Gamma(link = 'log')))

# using AIC, good case for location-scale models but
# I think for simplicity we stick with gamma models

AIC(result1[[1]], result[[1]])
AIC(result1[[6]], result[[6]])

# get approximate F-statistics and p-values
lapply(result, anova)

# estimated marginal means for interaction model
em1 <- function(x, div = 4:42){
  emmeans(x, ~ s(DIV, by = Cell.line), 
          type = 'response',
          at = list(DIV = div))
}



emm <- lapply(result, function(z) em1(z, div = 4:42))
emm


# plot time-series fixed effects



y_title <- expression(paste("Days ", italic("in vitro")))

plot.ts <- function(d, d0 = data_final_2, col.num = 1, ylab){
  pal <- c('#B8DE29FF',
           '#20A387FF',
           '#481567FF')
  ggplot(d %>% data.frame(), aes(x = DIV, y = response-1, group = Cell.line, colour = Cell.line)) +
    geom_point(data = d0, aes(y = d0[,col.num+4]),
               alpha = .2) +
    geom_ribbon(aes(ymin = lower.CL-1, ymax = upper.CL-1, fill = Cell.line, colour = NULL), alpha = .4) +
    geom_line(colour = 'black', size = 1) +
    facet_wrap(~Cell.line, 
               ncol = 1) +
    xlab(y_title) +scale_x_continuous(breaks = seq(0, 40, by = 5))+
    ylab(ylab) +
    theme_bw() +
    theme(axis.text.x=element_text( face="bold",size=15,color='black'),
          axis.text.y=element_text(size=14, face="bold",color='black',vjust=-0.5),
          axis.title.y=element_text(size=30, face="bold"),
          axis.title.x=element_text(size=30, face="bold"))+
    scale_fill_manual(values = pal, name = 'Cell line',labels = c('CLN3', 'CLN3-Cor'))+
    theme(legend.text=element_text(size=14),legend.title=element_text(size=16))
  
  ggsave(paste0("plot1_",col.num,".eps"), device=cairo_ps, width = 35, height = 23, units = "cm",fallback_resolution = 600)
  ggsave(paste0("plot1_",col.num,".png"), width = 15, height = 10, dpi = 300, units = "in", device='png')
  
}



plot2.ts <- function(d, d0 = data_final_2, col.num = 1, ylab){
  pal <- c('#808000',
           '#0D3B58')
  ggplot(d %>% data.frame(), aes(x = DIV, y = response-1, group = Cell.line, colour = Cell.line)) +
    geom_ribbon(aes(ymin = lower.CL-1, ymax = upper.CL-1, fill = Cell.line, colour = NULL), alpha = .4) +
    geom_line(colour = 'black', size = 1) +
    xlab(y_title) +
    ylab(ylab) +scale_x_continuous(breaks = seq(0, 40, by = 5))+
    theme_bw() +
    theme(axis.text.x=element_text( face="bold",size=15,color='black'),
          axis.text.y=element_text(size=14, face="bold",color='black',vjust=-0.5),
          axis.title.y=element_text(size=28, face="bold"),
          axis.title.x=element_text(size=30, face="bold"))+
    scale_fill_manual(values = pal, name = 'Cell line',labels = c('CLN3', 'CLN3-Cor'))+
    theme(legend.text=element_text(size=14),legend.title=element_text(size=16))+theme_classic()
  ggsave(paste0("plot2_",col.num,".eps"), device=cairo_ps, width = 35, height = 23, units = "cm",fallback_resolution = 600)
  ggsave(paste0("plot2_",col.num,".png"), width = 15, height = 10, dpi = 300, units = "in", device='png')
  
  
  
}






# y-axis labels
ylabs = c("Mean firing rate (Hz)",
          'Mean burst duration (s)',
          'Mean spike rate in burst (spikes/s)',
          '% spikes in burst',
          'Mean inter-burst interval (s)',
          'Mean network burst duration (s)',
          'Mean spike rate in network burst (s)',
          '% spikes in network burst',
          'Mean inter-network burst interval (s)',
          'Mean burst rate (bursts/min)',
          'Mean network burst rate(network bursts/min)',
          'Mean inter-spike interval (s)',
          '% active electrodes')

lapply(1:length(emm), function(z) plot.ts(emm[z],
                                          col.num = z,
                                          ylab = ylabs[z]))




lapply(1:length(emm), function(z) plot2.ts(emm[z],
                                           col.num = z,
                                           ylab = ylabs[z]))

# post-hoc contrasts
posthoc <- function(em){
  contrast(regrid(em), 
           'tukey',
           reverse = TRUE,
           interaction = list(DIV = 'identity',
                              Cell.line = 'tukey'))
}


ph <- lapply(emm, posthoc)
names(ph) <- ylabs
ph
sink()

# Principle components analysis using the method
# of singular-value decomposition
# scale each variable to mean = 0, SD = 1
pc <- prcomp(data_final_2[,cols], rank. = 1, scale. = TRUE)

# how much variance does the single PC explain?
summary(pc)


# fit a GAMM
data_final_2$PC1 <- pc$x

m1  <- mgcv::gam(PC1 ~s(DIV, by = Cell.line, bs = 'cs') + s(Well_Batch, bs = 're'),
                 data = data_final_2)
emm <- em1(m1)

pal <- c('#B8DE29FF',
         '#20A387FF')
g1 <- ggplot(emm %>% data.frame(), aes(x = DIV, y = emmean, group = Cell.line, colour = Cell.line)) +
  geom_point(data = data_final_2, aes(y = PC1),
             alpha = .2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Cell.line, colour = NULL), alpha = .4) +
  geom_line(colour = 'black', size = 1) +
  facet_wrap(~Cell.line, ncol = 1) +
  xlab('Days on MEA') +
  ylab('1st Principal Component') +
  theme_bw() +
  scale_colour_manual(values = pal, name = 'Cell line') +
  scale_fill_manual(values = pal, name = 'Cell line')

posthoc(emm)


# forest plot with 95% CI of the difference
# if the CI overlaps 0, the difference at that DIV 
# is not significant at p < .05 (two-sided)
g2 <- ggplot(posthoc(emm) %>% data.frame(),
             aes(x = DIV_identity, 
                 y = estimate,
                 colour = Cell.line_tukey)) +
  geom_linerange(aes(ymin = estimate-1.96*SE,
                     ymax = estimate+1.96*SE)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  coord_flip() +
  xlab('Days in-vitro') +
  ylab('Mean difference (z)') +
  scale_colour_manual(values = c('#B8DE29FF',
                                 '#20A387FF'),
                      name = '') 

cowplot::plot_grid(g1 + theme(legend.position = 'none'),
                   g2 + theme_bw() + theme(legend.position = 'bottom',
                                           legend.text = element_text(size = 6)),
                   labels = c('(a)', '(b)'),
                   label_x = -.02,
                   label_size = 10)
ggsave('pca_figure.png', dpi = 600, width = 6, height = 6)


library(stringr)
data<-read.csv("data_median2.csv", header=T, sep=",")
data2<-data %>% filter(print(across(everything(), ~ !str_detect(., "TOB220"))))
View(data2)
df<-data2
#scale
df[, -c(1)] <- scale(df[, -c(1)])


View(df)
#rownames(df) = make.names(df$Parameters, unique=TRUE)
df<-df[,-1]
View(df)
df_t<-t(df)
View(df_t)

#rename columns
colnames(df_t) <- c("1kb del_1", "1kb del_2", "1kb del_3", "1kb +_1", "1kb +_2", "1kb +_3")  
View(df_t)
df_matrix <- data.matrix(df_t)

View(df_t)



#replace NaN with zero

df_matrix[is.na(df_matrix)] = 0

View(df_matrix)

library(ggplot2)
library(ComplexHeatmap)
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
library(dendextend)
library(circlize)


b<-Heatmap(df_matrix,
           row_order = rownames(df_matrix),column_order =colnames(df_matrix),
           row_dend_width=unit(1.5, "cm"),
           column_names_gp = gpar(fontsize = 8),
           heatmap_legend_param = list(title = "Z score"),
           column_names_rot = 45,
           column_title = "",
           column_title_gp = gpar(fontsize = 20, fontface = "bold"),cluster_row_slices = FALSE,
           col=colorRamp2(c(-2,0,2),c("violet", "black", "springgreen")))
b

ggsave("heatmap_overview.eps", device=cairo_ps, width = 35, height = 26, units = "cm",fallback_resolution = 600)
ggsave("heatmap_overview.png",width = 15, height = 10, dpi = 300, units = "in", device='png')




