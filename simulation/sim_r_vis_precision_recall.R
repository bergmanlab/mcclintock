require(tools)
library(ggplot2)
library(plyr); library(dplyr)
library(tidyr)
library(Rmisc)
library(ggpubr)
library(reshape2)
# library(optparse)

# setup working directory
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript --vanilla sim_r_vis_precision_recall.R <mcc_dir> <sim_dir> <matrix_for_3_6_12_25_50_100>", call.=FALSE)
} else if (length(args) >= 3) {
      mcc_dir <- args[1]
      sim_dir <- args[2]
      table_list <- args[3:length(args)]
      # table_6x <- args[4]
      # table_12x <- args[5]
      # table_25x <- args[6]
      # table_50x <- args[7]
      # table_100x <- args[8]
}


# output dir
out_dir <- paste0(sim_dir, "/r_vis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## for interactive session
# mcc_dir <- "/home/jc33471/mcclintock/"
# sim_dir <- "/scratch/jc33471/mcc_paper/simulation_0912/"

##################### create table for plotting ##################### 

parse_table_list <- function(table_list){
  tables <- list()
  n_tables <- length(table_list)
  all_sim_table_covs <- data.frame()
  for (i in 1:n_tables){
      str <- unlist(strsplit(as.character(table_list[i]),"/")); str <- str[str!=""]
      cov <- as.numeric(str[length(str)-2])
      tables[[i]] <- read.csv(as.character(table_list[i]), stringsAsFactors = F) %>%
          select(contains("precision")|contains("recall")|contains("method")) %>%
          melt(id = "Method") %>%
          mutate(coverage = cov)
      all_sim_table_covs <- rbind(all_sim_table_covs, tables[[i]])
  }
  all_sim_table_covs <- all_sim_table_covs %>%
      separate(col = "variable", into = c("stat","window")) %>%
      mutate(Method = gsub("te-locate","te.locate",Method))
  all_sim_table_covs$window <- factor(all_sim_table_covs$window, levels = c("0","5","100","300","500"))
  return(all_sim_table_covs)
}
all_sim_table_covs <- parse_table_list(table_list)

##################### curve for precision and recall, all windows ##################### 

sim_precision_recall <- function(p_stat){
  
  get_plot_covs <- function(x){
    # modified title
    title <- c("ngs_te_mapper","RelocaTE","TEMP","RetroSeq","PoPoolationTE","TE-locate","ngs_te_mapper2","RelocaTE2","TEMP2","TEFLoN","PoPoolationTE2","TEBreak")
    names(title) <- c("ngs_te_mapper","relocate","temp","retroseq","popoolationte","te.locate","ngs_te_mapper2","relocate2","temp2","teflon","popoolationte2","tebreak")
    # set colors for lines
    setcolors <- c(
      "500" = "#cc4c02",
      "300" = "#ec7014",
      "100" = "#fe9929",
      "5" = "#fec44f",
      "0" = "#fee391"
    )
    labels <- c(
      "500" = "Within-500",
      "300" = "Within-300",
      "100" = "Within-100",
      "5" = "Within-5",
      "0" = "Exact" 
    )

    # create ggplot
    max_cov <- max(all_sim_table_covs$coverage)
    breaks_4 <- as.vector(quantile(0:max_cov))[-1]
    p <- all_sim_table_covs %>% filter(Method==x & stat==p_stat) %>% ggplot( aes(x=coverage, y=value, group=window, color=window)) +
      geom_line() +
      geom_point(aes(x=coverage,y=value)) +
      scale_color_manual(values = setcolors, labels = labels) +
      theme_bw() +
      scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_x_continuous(limits = c(0, max_cov), breaks = breaks_4) +
      labs(x = NULL, y = NULL, title = title[x]) +
      guides(color=guide_legend(title=NULL,nrow = 1)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14), text = element_text(size = 15))
    
    
    return(p)
  }
  
  p_ngs_te_mapper <- get_plot_covs("ngs_te_mapper")
  p_relocate <- get_plot_covs("relocate")
  p_temp <- get_plot_covs("temp")
  p_retroseq <- get_plot_covs("retroseq")
  p_popoolationte <- get_plot_covs("popoolationte")
  p_te.locate <- get_plot_covs("te.locate")
  p_ngs_te_mapper2 <- get_plot_covs("ngs_te_mapper2")
  p_relocate2 <- get_plot_covs("relocate2")
  p_temp2 <- get_plot_covs("temp2")
  p_teflon <- get_plot_covs("teflon")
  p_popoolationte2 <- get_plot_covs("popoolationte2")
  p_tebreak <- get_plot_covs("tebreak")
  
  # get legend in a separate file
#   legend <- ggpubr::get_legend(p_ngs_te_mapper)
#   ggpubr::as_ggplot(legend)
#   ggsave("~/Github/mcclintock2_paper/docs/manuscript/figures/legend.pdf", device = "pdf", width = 24, height = 0.5, units = "cm", dpi = 320)
  
  p_all <- ggpubr::ggarrange(p_ngs_te_mapper, p_relocate, p_temp, p_popoolationte, p_retroseq, p_te.locate,
                            p_ngs_te_mapper2, p_relocate2, p_temp2, p_popoolationte2, p_teflon, p_tebreak,
                            ncol = 6, nrow = 2, common.legend = TRUE, legend = "top")
  
  annotate_figure(p_all,
                  left = text_grob(p_stat, rot = 90, size = 25),
                  bottom = text_grob("Fold coverage", size = 15),
  )
  
  ggsave(paste0(out_dir,"/sim_",p_stat,"_allcovs.pdf"), device = "pdf", width = 13, height = 6, units = "in", dpi = 320)
}
sim_precision_recall("Precision")
sim_precision_recall("Recall")

##################### curve for precision and recall, 0 and 100 for main fig ##################### 

sim_precision_recall_main <- function(p_stat){
  
  get_plot_covs <- function(x){
    # modified title
    title <- c("ngs_te_mapper","RelocaTE","TEMP","RetroSeq","PoPoolationTE","TE-locate","ngs_te_mapper2","RelocaTE2","TEMP2","TEFLoN","PoPoolationTE2","TEBreak")
    names(title) <- c("ngs_te_mapper","relocate","temp","retroseq","popoolationte","te.locate","ngs_te_mapper2","relocate2","temp2","teflon","popoolationte2","tebreak")
    # set colors for lines
    setcolors <- c(
      "100" = "#E66100",
      "0" = "#5D3A9B"
    )
    labels <- c(
      "100" = "Within-100",
      "0" = "Exact" 
    )

    # create ggplot
    max_cov <- max(all_sim_table_covs$coverage)
    breaks_4 <- as.vector(quantile(0:max_cov))[-1]
    p <- all_sim_table_covs %>% filter(Method==x & stat==p_stat) %>% filter(window=="100" | window=="0") %>% ggplot( aes(x=coverage, y=value, group=window, color=window)) +
      geom_line() +
      geom_point(aes(x=coverage,y=value)) +
      scale_color_manual(values = setcolors, labels = labels) +
      theme_bw() +
      scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_x_continuous(limits = c(0, max_cov), breaks = breaks_4) +
      labs(x = NULL, y = NULL, title = title[x]) +
      guides(color=guide_legend(title=NULL,nrow = 1)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14), text = element_text(size = 15))
    
    
    return(p)
  }
  
  p_ngs_te_mapper <- get_plot_covs("ngs_te_mapper")
  p_relocate <- get_plot_covs("relocate")
  p_temp <- get_plot_covs("temp")
  p_retroseq <- get_plot_covs("retroseq")
  p_popoolationte <- get_plot_covs("popoolationte")
  p_te.locate <- get_plot_covs("te.locate")
  p_ngs_te_mapper2 <- get_plot_covs("ngs_te_mapper2")
  p_relocate2 <- get_plot_covs("relocate2")
  p_temp2 <- get_plot_covs("temp2")
  p_teflon <- get_plot_covs("teflon")
  p_popoolationte2 <- get_plot_covs("popoolationte2")
  p_tebreak <- get_plot_covs("tebreak")
  
  
  p_all <- ggpubr::ggarrange(p_ngs_te_mapper, p_relocate, p_temp, p_popoolationte, p_retroseq, p_te.locate,
                             p_ngs_te_mapper2, p_relocate2, p_temp2, p_popoolationte2, p_teflon, p_tebreak,
                             ncol = 6, nrow = 2, common.legend = TRUE, legend = "top")
  
  annotate_figure(p_all,
                  left = text_grob(p_stat, rot = 90, size = 25),
                  bottom = text_grob("Fold coverage", size = 15),
  )
  
  ggsave(paste0(out_dir,"/sim_",p_stat,"_mainfig.pdf"), device = "pdf", width = 13, height = 6, units = "in", dpi = 320)
}

sim_precision_recall_main("Precision")
sim_precision_recall_main("Recall")

