require(ggplot2)

my_theme_box <- theme_bw() + theme(
  text=element_text(size=20,  family="serif"),
  panel.grid=element_blank(), aspect.ratio=1,
  panel.background = element_blank())

my_theme_no_box <- my_theme_box + theme(
  panel.border = element_blank(),
  axis.line = element_line(colour = "black"))

my_theme <- my_theme_no_box

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

