
# circus plot --- 
## Figure4d. in vitro mutation rate for Label and Ctrl

library(tidyr)
library(dplyr)
library(circlize)
library(wesanderson)

setwd("./")

## 1. UG in label vs Ctrl

### A. creat the circle frame 

Convert_Circos_Coords = function(df){
  R_Loop = df %>%
    filter(Position < 577)
  Coding_Region = df %>%
    filter(Position > 576)%>%
    filter(Position < 16569)
  tail(R_Loop)
  head(Coding_Region)
  Circos_df = bind_rows(Coding_Region,R_Loop)
  Circos_df$Position = 1:length(Circos_df$Position)
  return(Circos_df)
}


Init_Circos_Mito = function(){
  pal <- wes_palette("Darjeeling1", 5, type = "discrete")
  circos.csv = read.csv("./Data/Circos_Table.csv")
  colnames(circos.csv) = c("Type", "Gene","Strand", "Start","End"  ) 
  mito_length = "16568"
  ticks = seq(0,17001,1000)
  df = data.frame(start = 0, end = 16568)
  rownames(df) = c("Mitochondria")
  
  circos.par("start.degree" = 90,"canvas.xlim" = c(-1.5,1.5))
  circos.par("gap.degree" = 0)
  circos.initialize(xlim = df)
  circos.par("track.height" = 0.3)
  circos.track(ylim = c(1,10))
  circos.axis(h = "top",
              major.at = ticks,
              labels.facing = "clockwise",
              major.tick.length = 2,
              labels.cex = 1)
  
  for(i in 1:length(circos.csv$Gene)){
    if(circos.csv$Strand[i] == '+'){
      if(circos.csv$Type[i] == "tRNA"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 5,col = pal[1])
        
      }
      if(circos.csv$Type[i] == "mRNA"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 5, col = pal[5])
        
      }
      if(circos.csv$Type[i] == "rRNA"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 5,col = "grey")
        
      }
      if(circos.csv$Type[i] == "NC"){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 10,ybottom = 1,col = "grey")
      }
    }
    if(circos.csv$Strand[i] == '-'){
      if(circos.csv$Type[i] == 'tRNA'){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 5,ybottom = 1,col = pal[1])
        
      }
      
      if(circos.csv$Type[i] == 'mRNA'){
        circos.rect(xleft = circos.csv$Start[i],xright = circos.csv$End[i],
                    ytop = 5,ybottom = 1,col = pal[5])
      }
    }
    if(circos.csv$End[i] < 8630){
      if(circos.csv$Type[i] == 'tRNA'){circos.text((circos.csv$Start[i]+circos.csv$End[i])/2, 7.5, '', facing = "clockwise", cex =0,font =2, niceFacing = T)
        }
      else{
      circos.text((circos.csv$Start[i]+circos.csv$End[i])/2, 7.5, circos.csv$Gene[i], facing = "clockwise", cex =0.8,font =2, niceFacing = T)
      }
      }else{
        if(circos.csv$Type[i] == 'tRNA'){circos.text((circos.csv$Start[i]+circos.csv$End[i])/2, 7.5, '', facing = "reverse.clockwise", cex =0,font =2)
        }
        else{
        
      circos.text((circos.csv$Start[i]+circos.csv$End[i])/2, 7.5, circos.csv$Gene[i], facing = "reverse.clockwise", cex = 0.8,font =2)
        }
    }
  }
  
  circos.text(4000, -23, 'Mutation rate', facing = "clockwise", cex = 1,font =1)
}


## 1. raw muttaion rate for TG set limit -----

### A. creat the circle frame 

circos.clear()
Init_Circos_Mito()
circos.par("track.height" = 0.4) 
circos.track(ylim = c(0,0.015))
circos.yaxis(side = c("left"),labels.niceFacing = TRUE,labels.cex = 1)


### B. input data and plot, raw mutation rate 

allinfo = read.table("./Data/vitro_mito_modrate_profiles.txt",header=TRUE,sep='\t')
head(allinfo)

#### lable
flt1 = allinfo %>% select(Position,Modified_rate)
colnames(flt1) = c("Position","mod_rate")
flt1 = Convert_Circos_Coords(flt1)
head(flt1)
median(flt1$mod_rate)
flt1$mod_rate[flt1$mod_rate >= 0.015] <- 0.015


#### Ctrl
flt2 = allinfo %>% select(Position,Untreated_rate)
colnames(flt2) = c("Position","mod_rate")
flt2 = Convert_Circos_Coords(flt2)
head(flt2)
median(flt2$mod_rate)
flt2$mod_rate[flt2$mod_rate >= 0.015] <- 0.015

### Plot Lines
circos.lines(x = 1:length(flt1$Position), y = flt1$mod_rate, lwd = 1, col = "#edae9a")
circos.lines(x = 1:length(flt2$Position), y = flt2$mod_rate, lwd = 1, col = "#e5e5e5")

