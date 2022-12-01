require(tools)
library(ggplot2)
library(plyr); library(dplyr)
library(tidyr)
library(Rmisc)
library(ggpubr)
library(reshape2)
pdf(NULL)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4 ) {
  stop("Usage: Rscript --vanilla sim_r_vis_for_cov.R <mcc_dir> <sim_dir> <cov> <TE_families>", call.=FALSE)
} else if (length(args) >= 4) {
    mcc_dir <- args[1]
    sim_dir <- args[2]
    cov <- args[3]
    te_list <- args[4:length(args)]
}

# output dir
out_dir <- paste0(sim_dir, "/r_vis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # path to simulation run
dir <- paste0(sim_dir, "/", cov)


# method list
methods <- c("ngs_te_mapper","relocate","temp","retroseq","popoolationte","te-locate","ngs_te_mapper2","relocate2","temp2","teflon","popoolationte2","tebreak")

##################### density plots ##################### 

# function to generate data files for each coverage
simulation_density_data <- function(te_list){
    # path to density .dist data
    output.location <- paste0(out_dir, "/density_data/",cov)
    dir.create(output.location, showWarnings = FALSE, recursive = TRUE)
    # simulation results location
    results.location <- paste0(dir, "/results/forward")
    reverse.results.location <- paste0(dir, "/results/reverse")
    
    # system(paste("cut -f1,2 ", mcc_dir,"/test/sacCer2.fasta.fai > ", output.location, "/sacCer2.lengths", sep=""))
    # prepare expected insertion site
    system(paste0("cat ", dir, "/data/forward/*.modref.bed | awk 'BEGIN{OFS=\"\t\"}{$2=$2+1; $3=$2; $5=0; $6=\"+\"; print $0}' | sort -k1,1 -k2,2n > ", output.location, "/singleinsertionlocationsTSS_for.bed"))
    system(paste0("cat ", dir, "/data/reverse/*.modref.bed | awk 'BEGIN{OFS=\"\t\"}{$2=$2+1; $3=$2; $5=0; $6=\"-\"; print $0}' | sort -k1,1 -k2,2n > ", output.location, "/singleinsertionlocationsTSS_rev.bed"))
    #system(paste0("cat ", output.location, "/singleinsertionlocationsTSS_for.bed ", output.location, "/singleinsertionlocationsTSS_rev.bed | sort -k1,1 -k2,2n > ", output.location, "/singleinsertionlocationsTSS.bed" ))
    for (te in te_list){
        system(paste("grep $'\t'", te, "$'\t' ", output.location, "/singleinsertionlocationsTSS_for.bed > ", output.location, "/singleinsertionlocationsTSS_for_", te, ".bed", sep=""))
        system(paste("grep $'\t'", te, "$'\t' ", output.location, "/singleinsertionlocationsTSS_rev.bed > ", output.location, "/singleinsertionlocationsTSS_rev_", te, ".bed", sep=""))
    }
    
    ## generate distance txt
    for(method in methods){
        all.bed <- vector()
        bedfile <- list.files(Sys.glob(paste0(results.location, "/run_*/*_1/results/", method)), full.names=T, pattern=paste0(sub("-","",method), "_nonredundant.bed$"))
        all.bed <- c(all.bed, bedfile)
        all.bed.rev <- vector()
        bedfile.rev <- list.files(Sys.glob(paste0(reverse.results.location, "/run_*/*_1/results/", method)), full.names=T, pattern=paste0(sub("-","",method), "_nonredundant.bed$"))
        all.bed.rev <- c(all.bed.rev, bedfile.rev)
        system(paste0("rm -f ", output.location, "/",method,".for.all.bed"))
        for(i in 1:length(all.bed))
        {
            system(paste("grep non- ", all.bed[i], " >> ", output.location, "/",method,".for.all.bed" ,sep=""))
        }
        system(paste0("rm -f ", output.location, "/",method,".rev.all.bed"))
        for(i in 1:length(all.bed.rev))
        {
            system(paste("grep non- ", all.bed.rev[i], " >> ", output.location, "/",method,".rev.all.bed" ,sep=""))
        }
        
        for(TE in te_list){
        # create file for each TE family
        system(paste("grep --no-filename ", TE, "\\| ", output.location, "/",method,".for.all.bed | grep non-reference | sort -k 1,1 > ", output.location, "/", TE, "_",method,".for.bed", sep=""))
        system(paste("grep --no-filename ", TE, "\\| ", output.location, "/",method,".rev.all.bed | grep non-reference | sort -k 1,1 > ", output.location, "/", TE, "_",method,".rev.bed", sep=""))
        # system(paste("~/miniconda/envs/mcc_plots/bin/bedtools genomecov -i ", output.location, "/", TE, "_",method,".bed -g ", output.location, "/sacCer2.lengths -bga -trackline | ~/miniconda/envs/mcc_plots/bin/wigToBigWig -clip stdin ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "_",method,".bw", sep=""))
        # find predictions with the window
        system(paste0("bedtools window -u -l 500 -r 505 -sw -a ",output.location,"/",TE,"_",method,".for.bed -b ",output.location,"/singleinsertionlocationsTSS_for_",TE,".bed | sort -k1,1 -k2,2n > ",output.location,"/",TE,"_",method,"_window.for.bed"))
        system(paste0("bedtools window -u -l 500 -r 505 -sw -a ",output.location,"/",TE,"_",method,".rev.bed -b ",output.location,"/singleinsertionlocationsTSS_rev_",TE,".bed | sort -k1,1 -k2,2n > ",output.location,"/",TE,"_",method,"_window.rev.bed"))
        # find closest predictions
        system(paste0("bedtools closest -D b -a ",output.location,"/",TE,"_",method,"_window.for.bed -b ",output.location,"/singleinsertionlocationsTSS_for_",TE,".bed > ",output.location,"/",TE,"_",method,"_dist.for.out"))
        system(paste0("bedtools closest -D b -a ",output.location,"/",TE,"_",method,"_window.rev.bed -b ",output.location,"/singleinsertionlocationsTSS_rev_",TE,".bed > ",output.location,"/",TE,"_",method,"_dist.rev.out"))
        # combine forward and reverse strands
        system(paste0("cat ",output.location,"/",TE,"_",method,"_dist.for.out ",output.location,"/",TE,"_",method,"_dist.rev.out > ", output.location, "/",TE,"_",method,"_dist.out"))
        # calculate and manipulate relative coordinates to expected TSS
        system(paste0("awk 'BEGIN{OFS=FS=\"\t\"}{chr=\"synthetic\";start=($2-$8+500);end=($3-$8+500); if(start >= 0 && end <= 1005) {print chr,start,end}}' ",output.location,"/",TE,"_",method,"_dist.out > ",output.location,"/",TE,"_",method,"_relative.bed"))
        # system(paste0("awk 'BEGIN{OFS=FS=\"\t\"}{chr=\"synthetic\";start=($2-$8+500);end=($3-$8+500);print chr,start,end}' ",output.location,"/",TE,"_",method,"_dist.rev.out > ",output.location,"/",TE,"_",method,"_relative.rev.bed"))
        # calculate bedgraph for plot
        system(paste0("echo \"synthetic\t1005\" > ",output.location,"/window.length"))
        num_syn <- system(paste0("cat ",output.location,"/singleinsertionlocationsTSS_for_",TE,".bed ", output.location,"/singleinsertionlocationsTSS_for_",TE,".bed | wc -l"), intern = T)
        system(paste0("bedtools genomecov -d -i ", output.location,"/",TE,"_",method,"_relative.bed -g ", output.location, "/window.length | awk 'BEGIN{OFS=\"\t\"}{$2=$2-500; $3=(",num_syn,"?$3/",num_syn,":0); print $0}'> ", output.location, "/",TE,"_",method,"_relative.bedgraph"))
        }
    }
}

## function to plot for method
simdensity_plots <- function(TE){
    # path to density .dist data
    output.location <- paste0(out_dir, "/density_data/",cov)
    dir.create(output.location, showWarnings = FALSE, recursive = TRUE)
    # methods title
    title2 <- c("ngs_te_mapper","RelocaTE","TEMP","RetroSeq","PoPoolationTE","TE-locate","TE-locate","ngs_te_mapper2","RelocaTE2","TEMP2","TEFLoN","PoPoolationTE2","TEBreak","Coverage")
    names(title2) <- c("ngs_te_mapper","relocate","temp","retroseq","popoolationte","te-locate","telocate","ngs_te_mapper2","relocate2","temp2","teflon","popoolationte2","tebreak","coverage")
    
    get_ty_plot <- function(method,TE){

        # tRNA profile plot
        dist_trna <- read.table(paste0(output.location, "/",TE,"_",method,"_relative.bedgraph"), stringsAsFactors = F)

        p1 <- ggplot(dist_trna, aes(x = V2,y = V3)) +
            xlim(c(-500,505)) +
            geom_line(aes(),color = "brown") +
            theme_bw() +
            #scale_fill_gradient(low = "white", high = "#0072B2") +
            theme(axis.text.x = element_text(), axis.ticks.x = element_line(),
                    axis.text.y = element_text(angle = 90, hjust = 0.5), axis.ticks.y = element_line(),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
            labs(x = NULL, y = NULL, title = title2[method])

        return(p1)
    }

    # plotSets <- list(); n <- 1
    # for (method in title2){
    #     plotSets[[n]] <- get_ty_plot(method,TE)
    #     n <- n + 1
    # }
    # generate plotset for all methods for $TE family
    plotSets <- lapply(methods,get_ty_plot,TE=TE)
    out_p <- ggpubr::ggarrange(plotlist = plotSets, nrow = 2, ncol = 6)
    out_p2 <- annotate_figure(out_p, left = text_grob(paste0(TE, " mean coverage"), rot = 90, size = 10) , bottom = text_grob("Position to insertion start (bp)", size = 10))
    ggsave(plot = out_p2,
        paste0(out_dir, "/simdensity_",cov,"x_",TE,".pdf"), device = "pdf", width = 20, height = 8, units = "in", dpi = 320)

}

simulation_density_data(te_list)

for (te in te_list){
    simdensity_plots(te)
}

##################### tsd plots ##################### 

simulation_tsd_len <- function(te_list){
    # settings for methods/title
    srmethods <- c("ngs_te_mapper","ngs_te_mapper2","relocate","relocate2","tebreak","temp")
    title2 <- c("ngs_te_mapper","RelocaTE","TEMP","RetroSeq","PoPoolationTE","TE-locate","TE-locate","ngs_te_mapper2","RelocaTE2","TEMP2","TEFLoN","PoPoolationTE2","TEBreak","Coverage")
    names(title2) <- c("ngs_te_mapper","relocate","temp","retroseq","popoolationte","te-locate","telocate","ngs_te_mapper2","relocate2","temp2","teflon","popoolationte2","tebreak","coverage")

    # path to density .dist data
    output.location <- paste0(out_dir, "/tsd_data/",cov)
    dir.create(output.location, showWarnings = FALSE, recursive = TRUE)
    results.location <- paste(dir, "/results/forward", sep="")
    reverse.results.location <- paste(dir, "/results/reverse", sep="")
    
    all.bed <- vector()
    for (method in srmethods){
        bedfile <- list.files(Sys.glob(paste0(results.location, "/run_*/*.modref_1/results/", method)), full.names=T, pattern=paste0(method, "_nonredundant.bed$"))
        bedfilerev <- list.files(Sys.glob(paste0(reverse.results.location, "/run_*/*.modref_1/results/", method)), full.names=T, pattern=paste0(method, "_nonredundant.bed$"))
        all.bed <- c(all.bed, bedfile, bedfilerev)
    }
    
    
    # For each bed file
    system(paste0("> ", output.location, "/all.bed"))
    for(i in 1:length(all.bed))
    {
        system(paste("grep non- ", all.bed[i], " >> ", output.location, "/all.bed" ,sep=""))
    }
    dataset <- read.table(paste(output.location, "/all.bed" ,sep=""), header=FALSE,col.names=c("Chromosome","Start","End", "TE Name", "Score", "Orientation"),sep="\t")
    
    dataset[grep("ngs_te_mapper",  dataset$TE.Name), "Method"] <- "ngs_te_mapper"
    dataset[grep("ngs_te_mapper2",  dataset$TE.Name), "Method"] <- "ngs_te_mapper2"
    dataset[grep("relocate",  dataset$TE.Name), "Method"] <- "relocate"
    dataset[grep("relocate2",  dataset$TE.Name), "Method"] <- "relocate2"
    dataset[grep("temp_sr",  dataset$TE.Name), "Method"] <- "temp"
    dataset[grep("temp_rp",  dataset$TE.Name), "Method"] <- "tempreadpair"
    dataset[grep("temp_nonab",  dataset$TE.Name), "Method"] <- "tempreference"
    dataset[grep("temp2_rp",  dataset$TE.Name), "Method"] <- "temp2readpair"
    dataset[grep("temp2_nonab",  dataset$TE.Name), "Method"] <- "temp2reference"
    dataset[grep("temp\\|sr",  dataset$TE.Name), "Method"] <- "temp"
    dataset[grep("temp\\|rp",  dataset$TE.Name), "Method"] <- "tempreadpair"
    dataset[grep("temp\\|nonab",  dataset$TE.Name), "Method"] <- "tempreference"
    dataset[grep("temp2\\|rp",  dataset$TE.Name), "Method"] <- "temp2readpair"
    dataset[grep("temp2\\|nonab",  dataset$TE.Name), "Method"] <- "temp2reference"
    dataset[grep("retroseq",  dataset$TE.Name), "Method"] <- "retroseq"
    dataset[grep("telocate",  dataset$TE.Name), "Method"] <- "telocate"
    dataset[grep("te-locate",  dataset$TE.Name), "Method"] <- "telocate"
    dataset[grep("popoolationte",  dataset$TE.Name), "Method"] <- "popoolationte"
    dataset[grep("popoolationte2",  dataset$TE.Name), "Method"] <- "popoolationte2"
    dataset[grep("teflon",  dataset$TE.Name), "Method"] <- "teflon"
    dataset[grep("tebreak",  dataset$TE.Name), "Method"] <- "tebreak"
    
    for (te in te_list){
        dataset[grep(paste0(te,"\\|"),  dataset$TE.Name), "TE"] <- te
    }

    dataset$tsd <- dataset$End - dataset$Start
    # new.TEs <- grep("non-reference",dataset$TE.Name)
    
    # Count the number of TEs per sample for each method
    # newTEs <- paste(dataset$Method[new.TEs],dataset$Sample[new.TEs])
    # new.TE.counts <- table(newTEs)
    get_ty_plot <- function(method,te){
        p1 <- dataset %>% filter(Method==method & TE==te) %>% ggplot(aes(x = tsd)) +
            xlim(c(0,15)) +
            geom_histogram(binwidth = 1, color="brown", fill="white") +
            theme_bw() +
            theme(axis.text.x = element_text(), axis.ticks.x = element_line(),
                    axis.text.y = element_text(angle = 90, hjust = 0.5), axis.ticks.y = element_line(),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
            scale_y_continuous(breaks = function(x) seq.int(ceiling(x[1]), floor(x[2]), by = ceiling((floor(x[2])-ceiling(x[1]))/4))) +
            labs(x = NULL, y = NULL, title = title2[method])
        return(p1)
    }

    for (TE in te_list){
        plotSets <- lapply(srmethods,get_ty_plot,te=TE)
        out_p <- ggpubr::ggarrange(plotlist = plotSets, nrow = 1)
        out_p2 <- annotate_figure(out_p, left = text_grob(paste0(TE, " count"), rot = 90, size = 10) , bottom = text_grob("Predicted TSD lengths (bp)", size = 10))
        ggsave(plot = out_p2,
            paste0(out_dir, "/simtsd_",cov,"x_",TE,".pdf"), device = "pdf", width = 20, height = 4, units = "in", dpi = 320)

    }
  
}
simulation_tsd_len(te_list)