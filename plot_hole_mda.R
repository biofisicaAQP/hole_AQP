# Required packages -------------------------------------------------------
# For this script "dplyr" is used for data processing."ggplot" is only used for final pot rendering.

library(ggplot2)
library(dplyr)
library(viridis)
library(hexbin)
library(reshape2)
library(DescTools)

# Load Function -----------------------------------------------------------
# Input dataset contains 4 columns: 'Time', 'Crap', 'Z coordinate' value and 'Radius' value.


load.dataset <- function(path){
  set <- read.table(path, sep="", header=FALSE)
  names(set) <- c("Time","Z", "Radius", "crap")
  set <- select(set,-c(crap))
  set <- set[-1,]
  print(paste0(format(object.size(set), units = "auto"), " - ", path))
  return(set)
}

all_chains <- data.frame(matrix(ncol = 0, nrow = 0))

for (i in c("cadenaA", "cadenaB", "cadenaC", "cadenaD")) #filenames (MDAnalysis output)
{
# Load dataset ------------------------------------------------------------
cadena <- i
ppath <-  paste(cadena,".txt",sep = "")
ds <- load.dataset(ppath)
ds.corr <- ds

j = "ProteinName"
secuencia <- seq(-30, 30, by=0.5) # Intervals
intervalos <- cut(ds.corr$Z, secuencia)
resumen.median <- tapply(ds.corr$Radius, intervalos, mean) # median or mean
res.tot <- data.frame(Zres=secuencia[1:(length(secuencia)-1)], 
                      Resumen.median=resumen.median,
                      Proteina=j)
rownames(res.tot) <- NULL

# Profile Plot ------------------------------------------------------------
# PlotLine <- ggplot(data=res.tot, aes(x=-Zres, y = Resumen.median, color=Proteina)) +
#   geom_line(linetype=1,size=1)+
#   xlab("Z coordinate (Å)") +
#   ylab("Radius (Å)")+
#   coord_cartesian(xlim = c(-30, 40), ylim = c(0, 4))+
#   scale_x_continuous(breaks=seq(-40,40,5))+
#   scale_y_continuous(breaks=seq(0,4,1))+
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   theme(legend.position="none")
# 
# ggsave(filename = paste(cadena,"_PoreProfile.png",sep = ""), plot = PlotLine, path = NULL,
#        scale = 1,
#        width = 25,
#        height = 15,
#        units = c("cm"),
#        dpi = 600)

# Heatmap plot ---------------------------------------------------------
  myColor = rev(RColorBrewer::brewer.pal(11, "Spectral"))
  myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans = "sqrt")
  
  
dat <- data.frame(x = ds$Z, y = ds$Radius)
p <- ggplot(dat, aes(x = x, y = y)) +
  myColor_scale_fill +
  stat_binhex(binwidth=c(0.25,0.025))+
  # geom_bin2d(binwidth=c(0.25,0.025)) +
  xlab("Z coordinate (Å)") +
  ylab("Radius (Å)")+
  coord_cartesian(xlim = c(-30, 20), ylim = c(0, 4))+
  scale_x_continuous(breaks=seq(-40,40,5))+
  scale_y_continuous(breaks=seq(0,4,1))+
  #scale_fill_viridis(option = "inferno")+
  #scale_fill_gradient(low="lightyellow", high="red1", trans = "sqrt")+ #mid="ivory"
  #geom_density2d(colour = "black")+
  #scale_fill_distiller(palette="RdOrYl")+
  geom_hline(yintercept=2, color = "red", size=0.5)+
  geom_hline(yintercept=1, color = "red", size=0.5)+
  geom_hline(yintercept=3, color = "red", size=0.5)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #theme(plot.title = element_text(size=30), axis.title=element_text(size=22, face="bold"))
 # scale_x_discrete(expand = c(0, 0))+
 # scale_y_continuous(expand = c(0,0))+
 
  
# theme(legend.position="none")

ggsave(filename = paste(cadena,"_heat.png",sep = ""), plot = p, path = NULL,
       scale = 1,
       width =27.2,
       height = 15,
       units = c("cm"),
       dpi = 200)

all_chains <- rbind(all_chains,ds)
}
dat <- data.frame(x = -all_chains$Z, y = all_chains$Radius)
p <- ggplot(dat, aes(x = x, y = y)) +
  stat_binhex(binwidth=c(0.25,0.025))+
  # geom_bin2d(binwidth=c(0.25,0.025)) +
  xlab("Z coordinate (Å)") +
  ylab("Radius (Å)")+
  coord_cartesian(xlim = c(-15, 45), ylim = c(0, 4))+
  scale_x_continuous(breaks=seq(-40,40,5))+
  scale_y_continuous(breaks=seq(0,4,1))+
  scale_fill_viridis(option = "inferno")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# theme(legend.position="none")

ggsave(filename = paste("all_chains_heat.png",sep = ""), plot = p, path = NULL,
       scale = 1,
       width =27.2,
       height = 15,
       units = c("cm"),
       dpi = 200)



j = "ProteinName"
secuencia <- seq(-30, 30, by=0.5) # Intervals
intervalos <- cut(all_chains$Z, secuencia)
resumen.median <- tapply(all_chains$Radius, intervalos, median) # median or mean
res.tot <- data.frame(Zres=secuencia[1:(length(secuencia)-1)], 
                      Resumen.median=resumen.median,
                      Proteina=j)
rownames(res.tot) <- NULL



# Profile Plot ------------------------------------------------------------
# PlotLine <- ggplot(data=res.tot, aes(x=-Zres, y = Resumen.median, color=Proteina)) +
#   geom_line(linetype=1,size=1)+
#   xlab("Z coordinate (Å)") +
#   ylab("Radius (Å)")+
#   coord_cartesian(xlim = c(-25, 25), ylim = c(0, 4))+
#   scale_x_continuous(breaks=seq(-25,25,5))+
#   scale_y_continuous(breaks=seq(0,4,1))+
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   theme(legend.position="none")
# 
# ggsave(filename = paste("all_chain_PoreProfile.png",sep = ""), plot = PlotLine, path = NULL,
#        scale = 1,
#        width = 25,
#        height = 15,
#        units = c("cm"),
#        dpi = 200)



# Min Radius -----------------------------------------------

# min.rads <- data.frame(unique(ds.corr$Time))
# names(min.rads) <- c("Frame")
# for (k in seq(10, 20)) {
#   z.low <- k
#   z.top <- z.low + 1
#   df <- data.frame()
#   for (i in unique(ds.corr$Time)) {
#     df <-
#       rbind(df, min(ds.corr$Radius[ds.corr$Z > z.low &
#                                               ds.corr$Z < z.top & ds.corr$Time == i]))
#   }
#   names(df) <- c(paste("Z",z.low,"-",z.top,sep = ""))
#   min.rads <- cbind(min.rads,df)
  # assign(paste("min.rad", z.low, z.top, sep = "."), df)
}
# data <- min.rads
# is.na(data)<-sapply(data, is.infinite)
# data[is.na(data)]<-NA
# d <- melt(data, id.vars = "Frame")
# PlotLine <- ggplot(d, aes(Frame,value, col=variable)) +
#   geom_line(aes(colour=),linetype=1,size=1) +
#   xlab("Frame") +
#   ylab("Radius (Å)")+
#   coord_cartesian(xlim = c(0, 2500), ylim = c(0, 5))+
#   scale_x_continuous(breaks=seq(0, 2500,500))+
#   scale_y_continuous(breaks=seq(0,5,1))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.position="none")
# PlotLine
# for (i in seq(2,12)){
# PlotSingleLine <- ggplot(min.rads, aes(x=min.rads$Frame, y=min.rads[,i])) +
#   geom_line(aes(colour=),linetype=1,size=1) +
#   xlab("Frame") +
#   ylab("Radius (Å)")+
#   coord_cartesian(xlim = c(0, 2500), ylim = c(0, 5))+
#   scale_x_continuous(breaks=seq(0, 2500,100))+
#   scale_y_continuous(breaks=seq(0,5,1))+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.position="none")
# print(PlotSingleLine)
# }
# 
# data <- min.rads
# is.na(data)<-sapply(data, is.infinite)
# data[is.na(data)]<-NA



quit()
