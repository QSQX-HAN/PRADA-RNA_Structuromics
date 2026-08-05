
# circus plot
## Extend Figure4r. in vitro reactivity

library(tidyr)
library(dplyr)
library(circlize)
library(wesanderson)
library(slider)

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
  
  circos.text(4000, -23, 'Reactivity', facing = "clockwise", cex = 1,font =1)
}


## 1. reactivity --- --- --- --- 

read.popavg = function(path,GU=TRUE){
  x = read.table(path,sep = "\t",header = T)
  if(GU){
    return(x)
  }
  else{
    x_AC = get_AC(x)
    return(x)
  }
}

get_UG = function(df){
  sub_df = df %>% filter(Sequence =="U" | Sequence =="G") 
  return(sub_df)
}

get_AC = function(df){
  sub_df = df %>%
    filter(Sequence == "A" | Sequence == "C") 
  return(sub_df)
}

DMS_React_Mean_Sliding_REAC = function(df1,window_size,UG=T,AC=F){
  if(UG && AC || !UG && !AC){
    print("Calculating Avg DMS reactivity on UGs and ACs")
    df1_Signal = df1
    Positions = df1$Position
  }
  if(UG){
    df1_Signal = get_UG(df1)
    Positions = df1_Signal$Position
  }
  if(AC){
    df1_Signal = get_AC(df1)
    Positions = df1_Signal$Position
  }
  Mus = df1_Signal$Norm_profile
  Mu_Windows = slide(Mus, ~.x, .after  = window_size)
  Pos_Windows = slide(Positions, ~.x, .after  = window_size)
  
  Mu_Mean_Vect = c()
  Position_Vector = c()
  length(Mu_Windows)
  for(i in 1:length(Mu_Windows)){
    
    Mu_Mean_Vect[i] = mean(Mu_Windows[[i]],na.rm=TRUE)
    Position_Vector[i] = median(Pos_Windows[[i]])
  }
  df_windows = tibble(Position = Position_Vector,
                      DMS_Mean = Mu_Mean_Vect)
  df_final = tibble(Position = df1$Position)
  df_final = left_join(df_final,df_windows)
  df_final = df_final %>%
    fill(DMS_Mean, .direction = "down")
  df_final
  return(df_final)
}


## Positive Strand Whole Genome

### Positive Strand Set Up
positive_coordinates = read.csv("./Data/MitoCoord_Pos.csv")
positive_coordinate_table = tibble(Position = c(),
                                   Type = c(),
                                   Gene = c())
for (i in 1:length(positive_coordinates$Gene)) {
  temp_tib  = tibble(Position = c(),
                     Type = c(),
                     Gene = c())
  Position = positive_coordinates$Start[i]:positive_coordinates$End[i]
  Type = rep(positive_coordinates$Type[i],length(Position))
  Gene = rep(positive_coordinates$Gene[i],length(Position))
  temp_tib = tibble(Position = Position,
                    Type = Type,
                    Gene = Gene)
  positive_coordinate_table = positive_coordinate_table %>% bind_rows(temp_tib)
}
rm(positive_coordinates)
rm(temp_tib)


pos_avg.list = list()
vitro_pos = read.popavg("./Data/vitro_mito_reactivity_profiles.txt")
#Add gene labels
vitro_pos = left_join(vitro_pos, positive_coordinate_table)
#Add to list
pos_avg.list[["vitro"]] = vitro_pos

#### 
vitro_Pos_react = DMS_React_Mean_Sliding_REAC(pos_avg.list$vitro,UG = T,AC = F,window_size = 80)

#### Convert Coordinates

vitro_Pos_react = Convert_Circos_Coords(vitro_Pos_react)

max(vitro_Pos_react$DMS_Mean)
min(vitro_Pos_react$DMS_Mean)

#### A. creat the circle frame 

circos.clear()
Init_Circos_Mito()
circos.par("track.height" = 0.4) 
circos.track(ylim = c(-0.2,3))
circos.yaxis(side = c("left"),labels.niceFacing = TRUE,labels.cex = 1)

#### Plot Lines
circos.lines(x = 1:length(vitro_Pos_react$Position), y = vitro_Pos_react$DMS_Mean, lwd = 2, col = "forestgreen")

