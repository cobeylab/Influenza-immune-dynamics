
 plot_themes  = 	theme_classic() +
  theme(axis.line = element_line(size=1)) +
  theme(axis.ticks = element_line(size=0.5)) +
  theme(axis.ticks.length = unit(-0.1,'cm')) +
  theme(axis.title.x=element_text(size=textSize)) +
  theme(axis.text.x=element_text(size= textSize, margin=margin(4,4,4,4,'pt'))) +
  theme(axis.title.y=element_text(size= textSize)) +
  theme(axis.text.y=element_text(size= textSize , margin=margin(4,4,4,4,'pt'))) +
  theme(plot.title=element_text(size=textSize)) +
  theme(legend.title=element_text(size=textSize)) +
  theme(legend.text=element_text(size=textSize)) +
  theme(legend.position ='right') +  theme(legend.direction='vertical') + theme(legend.box='horizontal')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.line = element_blank()) +
  theme(plot.margin = unit(c(.05, .05, .05, .05), "cm"))
