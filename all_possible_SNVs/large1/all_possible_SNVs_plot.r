
# Rscript --vanilla all_possible_SNVs_plot.r all_possible_SNVs

library(ggplot2)




args = commandArgs(trailingOnly=TRUE)
data = read.csv(args[1],header=T,sep="\t")



data_srz <- as.data.frame(table(data$classification))
level_order <- c('Synonymous', 'Missense', 'Nonsense', 'Start-loss')

p = ggplot(data_srz, aes(x =factor(Var1,level_order), y = Freq, fill = Var1)) + geom_bar(stat = "identity")  + geom_text(aes(label = Freq), vjust = -0.5, family="Decima Mono Pro", size=6) + theme_bw() + scale_fill_manual(values = c("Missense" = "lightblue","Nonsense"="red","Synonymous"="darkgreen", "Start-loss"="darkred"))




p = p + theme(axis.text.y=element_text(size=11,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
p = p + theme(axis.text.x=element_text(size=18,colour="black"),axis.title.x = element_blank())

p = p + theme(panel.grid.major=element_line(colour="white",size=0.5))
p = p + theme(panel.grid.minor=element_line(colour="white",size=0.5))

p = p + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))
p= p + theme(legend.text=element_text(colour="black",size=16))


p = p + theme(legend.position=c(0.8, 0.85),legend.title=element_text(colour="black",size=18))

p = p + guides(fill=guide_legend(title="Variant type"))

p = p + ylab('Variant Species') +labs(title="LARGE1 all possible SNVs")
p= p + theme(legend.background=element_rect(linetype="solid",colour="black"))


ggsave(filename=paste(basename(args[1]),".style1.png",sep = ""),dpi=300)
