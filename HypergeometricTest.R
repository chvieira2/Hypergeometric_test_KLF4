setwd("K:/Collaborations/Mauricio_Martins/KLF4_overlap_genes")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(limma)
library(ggrepel)
library(VennDiagram)
library(gplots)
library(data.table)




gene_table <- fread("summary comparison with single cell RNAseq.txt", sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = T)


gene_table[1, 3] <- 155





gene_table$p <- phyper(gene_table$overlap-1,
                        gene_table$cell.type.dataset,
                        23000-gene_table$cell.type.dataset,
                        gene_table$Klf4.dataset,
                        lower.tail = FALSE, log.p = FALSE)



gene_table$adjusted <-  p.adjust(gene_table$p, method="BH")







bar_plot <- ggplot(gene_table, aes(x=V1, y=enrichment)) +
  
  geom_bar(stat='identity', aes(fill = V1)) +
  
  annotate("text", x = "RGC", y = gene_table[gene_table$V1 == "RGC", "enrichment"],
           label = paste0("p.value <= ",round(gene_table[gene_table$V1 == "RGC", "p"], digits = 5)),
           hjust = 1,
           size = 7,
           fontface = "bold",
           color = rgb(0,0,0,1)) +
  
  coord_flip() + 
  ylab("Enrichment") +
  scale_x_discrete(limits=c("PRC","BP","HC","AC","MG","RGC")) +
  theme_bw() +
  theme(panel.border = element_blank(),
    legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(face="bold",
                                    size=25,
                                    hjust = 0.5,
                                    vjust = 1.5),
        axis.text.x  = element_text(face = "bold",
                                    color = "black",
                                    angle=0, 
                                    vjust=1,
                                    size=20),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(face = "bold",
                                    color = "black",
                                    angle=0, 
                                    vjust=.5,
                                    size=20))



png("barplot.png", width = 1000, height = 500, pointsize = 25)
grid.arrange(bar_plot, ncol=1)
dev.off()

rm(bar_plot)



names(gene_table) <- c("","cell.type", "cell.type.dataset", "Klf4.dataset", "overlap", "Genome", "enrichment", "p", "adjusted")

write.table(gene_table, file = "summary comparison with single cell RNAseq.txt", sep = "\t", na = "", quote = F)
