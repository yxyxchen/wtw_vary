library('ggplot2')
# plot theme 
saveTheme = theme(panel.background = element_rect(fill = "white")) + 
  theme(title =element_text(size = 20, face='bold'), 
        plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title=element_text(size= 20), axis.text = element_text(size=20), axis.line= element_line(color="black", size = 1)) +
  theme(strip.text = element_text(face="bold", size=15)) + 
  theme(legend.text=element_text(size= 15))

# plot theme 
displayTheme = theme(panel.background = element_rect(fill = "white",colour = "white"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+ 
  theme(axis.title=element_text(size=15), title =element_text(size=15, face='bold'), 
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size=15), axis.line= element_line(color="black", size = 0.5)) +
  theme(strip.text = element_text(face="bold", size=15)) + 
  theme(legend.text=element_text(size= 15))

myTheme = theme_linedraw() +
  theme(legend.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), text=element_text(face = "bold", size = 18),
        panel.border = element_rect(size = 1.5))

myThemeBig = theme_linedraw() +
  theme(legend.title = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), text=element_text(face = "bold", size = 23),
        panel.border = element_rect(size = 1.5))

symnum.args <- list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                    symbols = c("***", "**", "*", "ns"))
conditionColors =  c("#008837", "#7b3294")
paraColors = list("#e41a1c", "#377eb8", "#4daf4a", "#fe9929", "#ffd92f")