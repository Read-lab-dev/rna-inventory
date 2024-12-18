# Load library
library(VennDiagram)
library(dplyr)
DEG_RIOK2 <- read.csv("./bulk/DEseq2_GBM301_shRIOK2.csv")
DEG_IMP3 <- read.csv("./bulk/DEseq2_GBM301_sh77.csv")

select_RIOK2<- DEG_RIOK2 %>% filter(change=="up")
select_IMP3<- DEG_IMP3 %>% filter(change=="up")
select_RIOK2d<- DEG_RIOK2 %>% filter(change=="down")
select_IMP3d<- DEG_IMP3 %>% filter(change=="down")
# Generate 3 sets of 200 words
set1 <- select_RIOK2$X
set2 <- select_IMP3$X
set3 <- gene.sel
set4 <- select_RIOK2d$X
set5 <- select_IMP3d$X
# Prepare a palette of 3 colors with R colorbrewer:

library(RColorBrewer)
myCol <- brewer.pal(5, "Pastel2")[1:3]
myCol <- brewer.pal(5, "Pastel2")[5:3]
# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("shRIOK2-Up" , "shIMP3-Up" , "Unicox"),
  filename = 'venn_diagramm_up.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  col = "grey50",
  fill = myCol,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = .3,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.col = "black",
  cat.fontfamily = "sans",
  disable.logging	=T 
)

