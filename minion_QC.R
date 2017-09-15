#!/usr/bin/Rscript

args<-commandArgs(TRUE)

input.file = args[1]
output.dir = args[2]
q = as.numeric(as.character(args[3]))
cores = as.numeric(as.character(args[4]))

q_title = paste("Q>=", q, sep="")

# look at albacore's sequencing summary file
library(ggplot2)
library(viridis)
library(plyr)
library(reshape2)
library(yaml)
library(scales)
library(parallel)

# build the map for R9.5
p1 = data.frame(channel=33:64, row=rep(1:4, each=8), col=rep(1:8, 4))
p2 = data.frame(channel=481:512, row=rep(5:8, each=8), col=rep(1:8, 4))
p3 = data.frame(channel=417:448, row=rep(9:12, each=8), col=rep(1:8, 4))
p4 = data.frame(channel=353:384, row=rep(13:16, each=8), col=rep(1:8, 4))
p5 = data.frame(channel=289:320, row=rep(17:20, each=8), col=rep(1:8, 4))
p6 = data.frame(channel=225:256, row=rep(21:24, each=8), col=rep(1:8, 4))
p7 = data.frame(channel=161:192, row=rep(25:28, each=8), col=rep(1:8, 4))
p8 = data.frame(channel=97:128, row=rep(29:32, each=8), col=rep(1:8, 4))

q1 = data.frame(channel=1:32, row=rep(1:4, each=8), col=rep(16:9, 4))
q2 = data.frame(channel=449:480, row=rep(5:8, each=8), col=rep(16:9, 4))
q3 = data.frame(channel=385:416, row=rep(9:12, each=8), col=rep(16:9, 4))
q4 = data.frame(channel=321:352, row=rep(13:16, each=8), col=rep(16:9, 4))
q5 = data.frame(channel=257:288, row=rep(17:20, each=8), col=rep(16:9, 4))
q6 = data.frame(channel=193:224, row=rep(21:24, each=8), col=rep(16:9, 4))
q7 = data.frame(channel=129:160, row=rep(25:28, each=8), col=rep(16:9, 4))
q8 = data.frame(channel=65:96, row=rep(29:32, each=8), col=rep(16:9, 4))

map = rbind(p1, p2, p3, p4, p5, p6, p7, p8, q1, q2, q3, q4, q5, q6, q7, q8)

add_cols <- function(d, min.q){
    # repeat take a sequencing sumamry file (d),
    # return the same data frame with min.q and a cumulative.bases columns added
    d = subset(d, mean_qscore_template >= min.q)
    d = merge(d, map, by="channel")
    d = d[with(d, order(-sequence_length_template)), ] # sort by read length
    d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))
    d$hour = d$start_time %/% 3600
    
    # add the reads generated for each hour
    reads.per.hour = as.data.frame(table(d$hour))
    names(reads.per.hour) = c("hour", "reads_per_hour")
    reads.per.hour$hour = as.numeric(as.character(reads.per.hour$hour))
    d = merge(d, reads.per.hour, by = c("hour"))    
    return(d)
}


load_summary <- function(filepath, min.q){
    # load a sequencing summary and add some info
    # min.q is a vector of length 2 defining 2 levels of min.q to have
    d = read.delim(filepath)
    d$sequence_length_template = as.numeric(d$sequence_length_template)
    d$events_per_base = d$num_events_template/d$sequence_length_template
    
    flowcell = basename(dirname(filepath))
    
    d1 = add_cols(d, min.q[1])
    d1$Q_cutoff = "All reads"
    
    d2 = add_cols(d, min.q[2])
    d2$Q_cutoff = q_title
    
    d = rbind(d1, d2)
    d$flowcell = flowcell
    d$Q_cutoff = as.factor(d$Q_cutoff)
    return(d)
}

reads.gt <- function(d, len){
    return(length(which(d$sequence_length_template>=len)))
}

bases.gt <- function(d, len){
  reads = subset(d, sequence_length_template >= len)
  return(sum(as.numeric(reads$sequence_length_template)))
}

log10_minor_break = function (...){
    function(x) {
        minx         = floor(min(log10(x), na.rm=T))-1;
        maxx         = ceiling(max(log10(x), na.rm=T))+1;
        n_major      = maxx-minx+1;
        major_breaks = seq(minx, maxx, by=1)
        minor_breaks = 
            rep(log10(seq(1, 9, by=1)), times = n_major)+
            rep(major_breaks, each = 9)
        return(10^(minor_breaks))
    }
}

summary.stats <- function(d, Q_cutoff="All reads"){
    # Write summary stats for a single value of min.q
    rows = which(as.character(d$Q_cutoff)==Q_cutoff)
    d = d[rows,]
    
    total.bases = sum(as.numeric(d$sequence_length_template))
    N50.length = d$sequence_length_template[min(which(d$cumulative.bases > (total.bases/2)))]
    mean.length = round(mean(as.numeric(d$sequence_length_template)), digits = 1)
    median.length = round(median(as.numeric(d$sequence_length_template)), digits = 1)
    max.length = max(as.numeric(d$sequence_length_template))
    mean.q = round(mean(d$mean_qscore_template), digits = 1)
    median.q = round(median(d$mean_qscore_template), digits = 1)
    
    #calculate ultra-long reads and bases (max amount of data with N50>100KB)
    for(i in 1:nrow(d)){
        n = d$sequence_length_template[min(which(d$cumulative.bases>d$cumulative.bases[i]/2))]
        if(n < 100000){ break }
    }
    ultra.reads = as.integer(i-1)
    ultra.gigabases = sum(as.numeric(d$sequence_length_template[1:ultra.reads]))/1000000000
    
    reads = list(reads.gt(d, 20000), 
              reads.gt(d, 50000),
              reads.gt(d, 100000),
              reads.gt(d, 200000),
              reads.gt(d, 500000),
              reads.gt(d, 1000000),
              ultra.reads)
    names(reads) = c(">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m", "ultralong")

    bases = list(bases.gt(d, 20000)/1000000000, 
              bases.gt(d, 50000)/1000000000,
              bases.gt(d, 100000)/1000000000,
              bases.gt(d, 200000)/1000000000,
              bases.gt(d, 500000)/1000000000,
              bases.gt(d, 1000000)/1000000000,
              ultra.gigabases)
    names(bases) = c(">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m", "ultralong")
    
                
    return(list('total.gigabases' = total.bases/1000000000, 
                'N50.length' = N50.length, 
                'mean.length' = mean.length, 
                'median.length' = median.length,
                'max.length' = max.length,
                'mean.q' = mean.q,
                'median.q' = median.q,
                'reads' = reads,
                'gigabases' = bases
                ))
}

channel.summary <- function(d){
    a = ddply(d, .(channel), summarize, total.bases = sum(sequence_length_template), total.reads = sum(which(sequence_length_template>=0)), mean.read.length = mean(sequence_length_template), median.read.length = median(sequence_length_template))
    b = melt(a, id.vars = c("channel"))
    return(b)    
}

# supress warnings
options(warn=-1)



single.flowcell <- function(input.file, output.dir, q=8){

    print("Creating output directory")
    dir.create(output.dir)
    print(paste("Saving output to", output.dir))
    out.txt = file.path(output.dir, "summary.yaml")
    
    # write summaries
    print(paste("Loading and summarising input file, saving to:", out.txt))

    d = load_summary(input.file, min.q=c(-Inf, q))
    all.reads.summary = summary.stats(d, Q_cutoff = "All reads")
    q10.reads.summary = summary.stats(d, Q_cutoff = q_title)
    
    summary = list("input file" = input.file,
                   "All reads" = all.reads.summary,
                   cutoff = q10.reads.summary,
                   "notes" = 'ultralong reads refers to the largest set of reads with N50>100KB')
    
    names(summary)[3] = q_title
    
    write(as.yaml(summary), out.txt)
    
    # make plots
    print("Plotting length histogram")
    p1 = ggplot(d, aes(x = sequence_length_template)) + 
        geom_histogram(bins = 300) + 
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Read length") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "length_histogram.png"), width = 960/75, height = 960/75, plot = p1)

    print("Plotting mean Q score histogram")
    p2 = ggplot(d, aes(x = mean_qscore_template)) + 
        geom_histogram(bins = 300) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean Q score of read") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "q_histogram.png"), width = 960/75, height = 960/75, plot = p2)

    print("Plotting events per base histogram")
    p3 = ggplot(d, aes(x = events_per_base)) + 
        geom_histogram(bins = 300) + 
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean number of events per base called in read") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "epb_histogram.png"), width = 960/75, height = 960/75, plot = p3)
    
    print("Plotting flowcell overview")
    p5 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x=start_time/3600, y=sequence_length_template, colour = mean_qscore_template)) + 
        geom_point(size=1.5, alpha=0.35) + 
        scale_colour_viridis() + 
        labs(colour='Q')  + 
        scale_y_log10() + 
        facet_grid(row~col) +
        theme(panel.spacing = unit(0.5, "lines")) +
        xlab("Hours into run") +
        ylab("Read length") +
        theme(text = element_text(size = 40), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), legend.text=element_text(size=12))
    ggsave(filename = file.path(output.dir, "flowcell_overview.png"), width = 2500/75, height = 2400/75, plot = p5)

    print("Plotting flowcell yield summary")
    p6 = ggplot(d, aes(x=sequence_length_template, y=cumulative.bases, colour = Q_cutoff)) + 
        geom_line(size = 1) + 
        scale_x_continuous(breaks =c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000), limits = c(0, 100000)) +
        xlab("Minimum read length") +
        ylab("Total yield in bases") +
        scale_colour_discrete(guide = guide_legend(title = "Reads")) +
        theme(text = element_text(size = 15))
    ggsave(filename = file.path(output.dir, "yield_summary.png"), width = 960/75, height = 960/75, plot = p6)
    
    
    print("Plotting sequence length over time")
    e = subset(d, Q_cutoff=="All reads")
    e$Q = paste(">=", q, sep="")
    e$Q[which(e$mean_qscore_template<q)] = paste("<", q, sep="")
    p7 = ggplot(e, aes(x=start_time/3600, y=sequence_length_template, colour = Q, group = Q)) + 
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean read length") + 
        ylim(0, NA)
    ggsave(filename = file.path(output.dir, "length_by_hour.png"), width = 960/75, height = 480/75, plot = p7)
    
    
    print("Plotting Q score over time")
    p8 = ggplot(e, aes(x=start_time/3600, y=mean_qscore_template, colour = Q, group = Q)) + 
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean Q score") + 
        ylim(0, NA)
    ggsave(filename = file.path(output.dir, "q_by_hour.png"), width = 960/75, height = 480/75, plot = p8)
    

    print("Plotting number of reads over time")
    f = d[c("hour", "reads_per_hour", "Q_cutoff")]
    f = f[!duplicated(f),]
    g = subset(f, Q_cutoff=="All reads")
    h = subset(f, Q_cutoff==q_title)
    max = max(f$hour)
    all = 0:max
    add.g = all[which(all %in% g$hour == FALSE)]
    if(length(add.g)>0){
        add.g = data.frame(hour = add.g, reads_per_hour = 0, Q_cutoff = "All reads")
        g = rbind(g, add.g)
    }
    
    add.h = all[which(all %in% h$hour == FALSE)]
    if(length(add.h)>0){
        add.h = data.frame(hour = add.h, reads_per_hour = 0, Q_cutoff = q_title)
        h = rbind(h, add.h)
    }
    
    i = rbind(g, h)
    i$Q_cutoff = as.character(i$Q_cutoff)
    i$Q_cutoff[which(i$Q_cutoff==q_title)] = paste("Q>=", q, sep="")
    p9 = ggplot(i, aes(x=hour, y=reads_per_hour, colour = Q_cutoff, group = Q_cutoff)) + 
        geom_point() +
        geom_line() +
        xlab("Hours into run") + 
        ylab("Number of reads per hour") + 
        ylim(0, NA) + 
        scale_color_discrete(guide = guide_legend(title = "Reads"))
    ggsave(filename = file.path(output.dir, "reads_per_hour.png"), width = 960/75, height = 480/75, plot = p9)
    
        
    print("Plotting read length vs. q score scatterplot")
    p10 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + 
        geom_point(alpha=0.05, size = 0.4) + 
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        scale_colour_viridis(trans = "log", labels = scientific) + 
        labs(colour='Events per base\n(log scale)\n')  + 
        theme(text = element_text(size = 15)) +
        xlab("Read length") +
        ylab("Mean Q score of read")
    ggsave(filename = file.path(output.dir, "length_vs_q.png"), width = 960/75, height = 960/75, plot = p10)
    
    print("Plotting flowcell channels summary histograms")
    c = channel.summary(subset(d, Q_cutoff=="All reads"))
    c10 = channel.summary(subset(d, Q_cutoff==q_title))
    c$Q_cutoff = "All reads"
    c10$Q_cutoff = q_title
    cc = rbind(c, c10)
    cc$variable = as.character(cc$variable)
    cc$variable[which(cc$variable=="total.bases")] = "Number of bases per channel"
    cc$variable[which(cc$variable=="total.reads")] = "Number of reads per channel"
    cc$variable[which(cc$variable=="mean.read.length")] = "Mean read length per channel"
    cc$variable[which(cc$variable=="median.read.length")] = "Median read length per channel"
    
    p11 = ggplot(cc, aes(x = value)) + geom_histogram(bins = 30) + 
        facet_grid(Q_cutoff~variable, scales = "free_x") + 
        theme(text = element_text(size = 20))
    ggsave(filename = file.path(output.dir, "channel_summary.png"), width = 2400/75, height = 960/75, plot = p11) 
    
    # add the flowcell to d
    flowcell = basename(dirname(input.file))
    
    return(d)
}

multi.flowcell = function(input.file, output.base, q){
    
    print(paste("Creating folder", output.base))
    dir.create(output.base)
    flowcell = basename(dirname(input.file))
    dir.create(file.path(output.base, flowcell))
    output.dir = file.path(output.base, flowcell, "minionQC")
    d = single.flowcell(input.file, output.dir, q)
    
}


multi.plots = function(dm, output.dir){
    
    # make plots
    print("")
    print("**** Making combined plots ****")
    print("Plotting multi length distributions")
    p1 = ggplot(dm, aes(x = sequence_length_template)) + 
        geom_line(stat="density", aes(colour = flowcell, y = ..count..)) +
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Read length") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "length_distributions.png"), width = 960/75, height = 960/75, plot = p1)
    
    print("Plotting mean Q score distributions")
    p2 = ggplot(dm, aes(x = mean_qscore_template)) + 
        geom_line(stat="density", aes(colour = flowcell, y = ..count..)) +
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean Q score of read") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "q_distributions.png"), width = 960/75, height = 960/75, plot = p2)
    
    print("Plotting events per base distributions")
    p3 = ggplot(dm, aes(x = events_per_base)) + 
        geom_line(stat="density", aes(colour = flowcell, y = ..count..)) +
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        theme(text = element_text(size = 15)) +
        xlab("Mean number of events per base called in read") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "epb_distributions.png"), width = 960/75, height = 960/75, plot = p3)
    

    print("Plotting flowcell yield summary")
    p6 = ggplot(dm, aes(x=sequence_length_template, y=cumulative.bases, colour = flowcell)) + 
        geom_line(size = 1) + 
        scale_x_continuous(breaks =c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000), limits = c(0, 100000)) +
        xlab("Minimum read length") +
        ylab("Total yield in bases") +
        theme(text = element_text(size = 15)) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y")
    ggsave(filename = file.path(output.dir, "yield_summary.png"), width = 960/75, height = 960/75, plot = p6)
    

    print("Plotting sequence length over time")
    e = subset(dm, Q_cutoff=="All reads")
    e$Q = paste(">=", q, sep="")
    e$Q[which(e$mean_qscore_template<q)] = paste("<", q, sep="")
    p7 = ggplot(e, aes(x=start_time/3600, y=sequence_length_template, colour = flowcell)) + 
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean read length") + 
        ylim(0, NA) +
        facet_wrap(~Q, ncol = 1, scales = "free_y")
    ggsave(filename = file.path(output.dir, "length_by_hour.png"), width = 960/75, height = 480/75, plot = p7)
    

    print("Plotting Q score over time")
    p8 = ggplot(e, aes(x=start_time/3600, y=mean_qscore_template, colour = flowcell)) + 
        geom_smooth() + 
        xlab("Hours into run") + 
        ylab("Mean Q score") + 
        facet_wrap(~Q, ncol = 1, scales = "free_y")
    ggsave(filename = file.path(output.dir, "q_by_hour.png"), width = 960/75, height = 480/75, plot = p8)
    
        
    print("Plotting read length vs. q score scatterplot")
    point.size = 0.02 / length(unique(dm$flowcell))
    point.alpha = 0.04 / (length(unique(dm$flowcell)) * 0.5)
    p9 = ggplot(subset(dm, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + 
        geom_point(alpha=point.alpha, size = point.size) + 
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        scale_colour_viridis(trans = "log", labels = scientific) + 
        labs(colour='Events per base\n(log scale)\n')  + 
        theme(text = element_text(size = 15)) +
        xlab("Read length") +
        ylab("Mean Q score of read") + 
        facet_wrap(~flowcell)
    ggsave(filename = file.path(output.dir, "length_vs_q.png"), width = 960/75, height = 960/75, plot = p9)
    

}



if(file_test("-f", input.file)==TRUE){
    # if it's an existing file (not a folder) just run one analysis
    d = single.flowcell(input.file, output.dir, q)
}else if(file_test("-d", input.file)==TRUE){
    # it's a directory, recursively analyse all sequencing_summary.txt files
  
    # get a list of all sequencing_summary.txt files, recursively
    summaries = list.files(path = input.file, pattern = "sequencing_summary.txt", recursive = TRUE, full.names = TRUE)
      
    # analyse each one and keep the returns in a list
    print("")
    print("**** Analysing each input file from this list ****")
    print(summaries)
    results = mclapply(summaries, multi.flowcell, output.dir, q, mc.cores = cores)

    # rbind that list
    dm = do.call("rbind", results)

    multi.plots(dm, output.dir)
    
}else{
    #WTF
    warning("Could find a sequencing summary file in your input. The input must be either a sequencing_summary.txt file, or a directory containing one or more such files")
}
