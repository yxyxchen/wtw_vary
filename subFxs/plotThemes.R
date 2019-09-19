library('ggplot2')
# 
themeColor = "#3876AE"

# annotations for p values 
symnum.args <- list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                    symbols = c("***", "**", "*", "ns"))
# plot theme
myTheme = theme(legend.position = 'none', panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), text=element_text(face = "bold", size = 18),
                axis.line = element_line(size = 1),
                panel.background = element_rect(fill = "white"))
