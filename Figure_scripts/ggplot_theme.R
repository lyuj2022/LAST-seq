# ggplot theme
theme_man <- 
  theme_classic() +
  theme(
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(colour='black', size=7),
    legend.title=element_blank(),
    legend.key=element_blank(),
    axis.line.x = element_line(colour = "black", size=0.25),
    axis.line.y = element_line(colour = "black", size=0.25),
    axis.ticks.x = element_line(colour = "black", size = 0.25),
    axis.ticks.y = element_line(colour = "black", size = 0.25),
    strip.background=element_blank(),
    strip.text=element_text(size=7))

#rename x and y axis
xlab("UMI counts")
  ylab("cell counts")
  
#change x and y axis number formar with scales
library(scales)
  
  number_format(accuracy = NULL, scale = 1, prefix = "", suffix = "",
                big.mark = " ", decimal.mark = ".", trim = TRUE, ...)
  
  number(x, accuracy = NULL, scale = 1, prefix = "", suffix = "",
         big.mark = " ", decimal.mark = ".", trim = TRUE, ...)
  
  comma_format(accuracy = NULL, scale = 1, prefix = "", suffix = "",
               big.mark = ",", decimal.mark = ".", trim = TRUE, digits, ...)
  
  comma(x, accuracy = NULL, scale = 1, prefix = "", suffix = "",
        big.mark = ",", decimal.mark = ".", trim = TRUE, digits, ...)
  
  percent_format(accuracy = NULL, scale = 100, prefix = "",
                 suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE,
                 ...)
  
  percent(x, accuracy = NULL, scale = 100, prefix = "",
          suffix = "%", big.mark = " ", decimal.mark = ".", trim = TRUE,
          ...)
  
  e.g.
  
  scale_x_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ","))

  #replace scientific format with 10^..
  scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  }
  # the parse()函数能将字符串转换为表达式expression
  # The gsub() function in R is used for replacement operations. The functions takes the input and substitutes it against the specified values.
  # %*% = x
#############lenged
theme(legend.position = c(0.8, 0.2),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1, "cm"),
      legend.text = element_text(size=10))
#remove legend
theme(legend.position='none')
#rename legend
scale_fill_discrete(labels=c("LASTseq","SMARTseq"))
scale_color_discrete(labels=c("LASTseq","SMARTseq"))

#######add text
annotation_logticks(sides = "bl", size=0.15, short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
annotate(geom="text",x=10, y=1000, label="R^2=0.86",size=2)



#########save 
#width=45mm is one fourth of A4
ggsave(filename = "fittingeg2.pdf",plot = p,width=45,height = 45,units = "mm", dpi = 300)