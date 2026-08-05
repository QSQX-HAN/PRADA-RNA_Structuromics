
# circus plot --- 
## Figure4c. in vitro coverage

library(tidyr)
library(dplyr)
library(circlize)
library(wesanderson)

setwd("./")

## 1. coverage --- ---

### functions used 

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
  
  circos.text(4000, -23, 'Read depth', facing = "clockwise", cex = 1,font =1)
}

### A. input data and plot, coverage 

allinfo = read.table("./Data/vitro_circus_coverage.txt",header=TRUE,sep='\t')
head(allinfo)

flt1 = allinfo %>% select(Position,Coverage)
flt1$Position <- unlist(flt1$Position)
Cov_Circos = Convert_Circos_Coords(flt1)
head(flt1)
max(flt1$Coverage)

max_positive_coverage = max(Cov_Circos$Coverage)
max_positive_coverage = max_positive_coverage[1]

### B. start the cirus plot

circos.clear()
Init_Circos_Mito()
circos.par("track.height" = 0.4) 

gap = 0
circos.track(ylim = c(0,(max_positive_coverage+gap)))
circos.yaxis(
  side = c("left"),
  labels.niceFacing = TRUE,
  labels.cex = 0.5,
  at = c(-50000,-150000,-250000,-350000,
         50000,250000,450000,650000,850000
  )
)

threshold = 500
for (i in 1:length(Cov_Circos$Position)) {
  
  if(Cov_Circos$Coverage[i] > threshold){
    circos.rect(xleft = Cov_Circos$Position[i] ,xright = (Cov_Circos$Position[i] + 1),
                ybottom = gap, ytop =Cov_Circos$Coverage[i]+gap,col = "#8E28A4", border ="#8E28A4")
  }else{
    circos.rect(xleft = Cov_Circos$Position[i] ,xright = (Cov_Circos$Position[i] + 1),
                ybottom = gap, ytop =Cov_Circos$Coverage[i]+gap,col = "#C0D5F9", border ="#C0D5F9")
  }
}
