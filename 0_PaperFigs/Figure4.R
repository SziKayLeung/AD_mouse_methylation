
library("ggplot2")
library("cowplot")
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))
zenDir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/"


# manhattan plots
load(file = paste0(zenDir, "/2_annotated/manhattan_cumu_ls.RData"))
pman1 <- plot_grid(plot_manhattan_3(cumu_ls$rTg4510_ECX_interaction,"ECXRRBS","Tg4510","no"),
          plot_manhattan_3(cumu_ls$rTg4510_ECX_pathology,"ECXRRBS_upwards","Tg4510","no") + scale_y_reverse(lim=c(320,0)),
          ncol=1)
pman2 <- plot_grid(plot_manhattan_3(cumu_ls$J20_ECX_interaction,"ECXRRBS","J20","no"),
          plot_manhattan_3(cumu_ls$J20_ECX_pathology,"ECXRRBS_upwards","J20","no") + scale_y_reverse(lim=c(17,0)), ncol = 1)


# epigenetic clocks
load(file = paste0(zenDir, "/5_epigeneticClock/rTg4510Clock.RData"))
load(file = paste0(zenDir, "/5_epigeneticClock/J20Clock.RData"))

plot_grid(
  plot_clock(rTg4510Clocks$ECX,"ECX",model="rTg4510"),
  plot_clock(rTg4510Clocks$HIP,"HIP",model="rTg4510"),
  plot_clock(J20Clocks$ECX,"ECX",model="J20"),
  plot_clock(J20Clocks$HIP,"HIP",model="J20") , labels = c("A","B","C","D")
)

plot_grid(
  plot_clock(rTg4510Clocks$ECX,"ECX",model="rTg4510", boxplot = FALSE),
  plot_clock(rTg4510Clocks$HIP,"HIP",model="rTg4510", boxplot = FALSE),
  plot_clock(J20Clocks$ECX,"ECX",model="J20", boxplot = FALSE),
  plot_clock(J20Clocks$HIP,"HIP",model="J20", boxplot = FALSE),
  labels = c("A","B","C","D")
)

plot_grid(
  clock_acceleration(rTg4510Clocks$ECX, "rTg4510", "ECX"),
  clock_acceleration(rTg4510Clocks$HIP, "rTg4510", "HIP"),
  clock_acceleration(J20Clocks$ECX, "J20", "ECX"),
  clock_acceleration(J20Clocks$HIP, "J20", "HIP"),
  labels = c("A","B","C","D")
)

# 2_simple_correlation.R
As3mtimage <- lapply(As3mtimage, function(x) x + scale_colour_manual(values = c("#00AEC9","black")) + 
                       theme(strip.background=element_rect(colour="white", fill="white"),
                             legend.position = "None"))

As3mtimage1 <- lapply(As3mtimage, function(x) x + labs(x = NULL, y = NULL) + theme(legend.position = "None"))
As3mtimageListPlots <- plot_grid(As3mtimage1[[1]], As3mtimage1[[2]],
                                 As3mtimage1[[4]] + scale_x_continuous(breaks=seq(70, 100, 10)), 
                                 As3mtimage1[[8]], nrow=2)


y.grob <- textGrob("RNA-Seq normalized counts", gp=gpar(fontsize=12), rot=90)
x.grob <- textGrob("Methylation (%)", gp=gpar(fontsize=12))
pAs3mt <- grid.arrange(arrangeGrob(As3mtimageListPlots, left = y.grob, bottom = x.grob), nrow = 1)

pClocks <- plot_grid(pClock1, pClock2, labels = c("i","ii"), scale = 0.95)
allP <- plot_grid(
  plot_grid(pman1,pman2, labels = c("A","B"), scale = 0.95),
  plot_grid(pAs3mt, labels = c("C"), scale = 0.95),
  plot_grid(NULL, pClocks, labels = c("D","E"), scale = 0.95, nrow=1), ncol = 1, rel_heights = c(0.4,0.4,0.2), scale = 0.95
)

pdf(paste0("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/Figures/MainFigures4.pdf"), width = 10, height = 15)
allP
dev.off()