#!/usr/bin/env

args<-commandArgs(TRUE)

input.file = args[1]
output.dir = args[2]

# look at albacore's sequencing summary file
library(ggplot2)
library(viridis)
library(plyr)
library(reshape2)
library(yaml)
library(ggjoy)
library(scales)
library(data.table)

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
    setDT(d)[, reads_per_hour := sum(which(sequence_length_template>=1)), by = hour] 
    return(d)
}


load_summary <- function(filepath, min.q=c(-Inf, 10)){
    # load a sequencing summary and add some info
    # min.q is a vector of length 2 defining 2 levels of min.q to have
    d = read.delim(filepath)
    d$sequence_length_template = as.numeric(d$sequence_length_template)
    d$events_per_base = d$num_events_template/d$sequence_length_template
    
    flowcell = basename(dirname(filepath))
    
    d1 = add_cols(d, min.q[1])
    d1$Q_cutoff = "All reads"
    
    d2 = add_cols(d, min.q[2])
    d2$Q_cutoff = "Reads with mean Q score > 10"
    
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


print("Creating output directory")
dir.create(output.dir)
out.txt = file.path(output.dir, "summary.yaml")

# write summaries
print("Loading and summarising input file")
d = load_summary(input.file)
all.reads.summary = summary.stats(d, Q_cutoff = "All reads")
q10.reads.summary = summary.stats(d, Q_cutoff = "Reads with mean Q score > 10")

summary = list("input file" = input.file,
               "All reads" = all.reads.summary,
               "Reads with Q>10" = q10.reads.summary,
               "notes" = 'ultralong reads refers to the largest set of reads with N50>100KB')

write(as.yaml(summary), out.txt)

# make plots
print("Plotting length histogram")
png(filename = file.path(output.dir, "length_histogram.png"), width = 960, height = 960)
ggplot(d, aes(x = sequence_length_template)) + geom_histogram(bins = 500) + scale_x_log10() + facet_wrap(~Q_cutoff, ncol = 1) + theme(text = element_text(size = 15))
dev.off()

print("Plotting mean Q score histogram")
png(filename = file.path(output.dir, "q_histogram.png"), width = 960, height = 960)
ggplot(d, aes(x = mean_qscore_template)) + geom_histogram(bins = 500) + facet_wrap(~Q_cutoff, ncol = 1) + theme(text = element_text(size = 15))
dev.off()

print("Plotting events per base histogram")
png(filename = file.path(output.dir, "epb_histogram.png"), width = 960, height = 960)
ggplot(d, aes(x = events_per_base)) + geom_histogram(bins = 500) + scale_x_log10() + facet_wrap(~Q_cutoff, ncol = 1) + theme(text = element_text(size = 15))
dev.off()

print("Plotting flowcell channels events per base plot")
png(filename = file.path(output.dir, "flowcell_channels_epb.png"), width = 2400, height = 2400)
ggplot(subset(d, Q_cutoff=="All reads"), aes(x=start_time/3600, y=events_per_base, colour = mean_qscore_template)) + 
    geom_point(size=1.5, alpha=0.35) + 
    scale_colour_viridis() + 
    scale_y_log10(limits=c(1e+0, 1e+01)) + 
    facet_grid(row~col) +
    theme(panel.spacing = unit(0.5, "lines")) +
    xlab("Hours into run") +
    theme(legend.position="none") +
    theme(text = element_text(size = 40), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12))
dev.off()


print("Plotting flowcell yield summary")
png(filename = file.path(output.dir, "yield_summary.png"), width = 960, height = 960)
ggplot(d, aes(x=sequence_length_template, y=cumulative.bases, colour = Q_cutoff)) + 
    geom_line(size = 1) + 
    scale_x_continuous(breaks =c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000), limits = c(0, 100000)) +
    xlab("Minimum read length") +
    ylab("Total yield in bases") +
    theme(text = element_text(size = 15))
dev.off()


print("Plotting sequence length over time")
png(filename = file.path(output.dir, "length_by_hour.png"), width = 960, height = 960)
upper.limit = signif(quantile(subset(d, Q_cutoff=="All reads")[,sequence_length_template], probs = 0.95), 1)
ggplot(subset(d, sequence_length_template >= 1), aes(x = sequence_length_template, y = hour, group = hour, fill = reads_per_hour)) + 
    geom_joy(scale = 2, alpha = 0.8) + 
    scale_x_continuous(expand = c(0.01, 0), limits = c(0, upper.limit), labels = scientific) + 
    scale_y_reverse(expand = c(0.01, 0)) + 
    facet_wrap(~Q_cutoff, ncol = 1) +
    theme(text = element_text(size = 15)) +
    scale_fill_viridis()
dev.off()

print("Plotting Q score over time")
png(filename = file.path(output.dir, "q_by_hour.png"), width = 960, height = 960)
ggplot(d, aes(x = mean_qscore_template, y = hour, group = hour, fill = reads_per_hour)) + 
    geom_joy(scale = 3, alpha = 0.8) + 
    scale_x_continuous(expand = c(0.01, 0)) + 
    scale_y_reverse(expand = c(0.01, 0)) + 
    facet_wrap(~Q_cutoff, ncol = 1) +
    theme(text = element_text(size = 15)) +
    scale_fill_viridis()
dev.off()

print("Plotting read length vs. q score scatterplot")
d$events_per_base[which(d$events_per_base>10)] = 10
png(filename = file.path(output.dir, "length_vs_q.png"), width = 960, height = 960)
ggplot(subset(d, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + geom_point(alpha=0.05, size = 0.4) + scale_x_log10(breaks=c(1e+01, 1e+02, 1e+03,1e+04,1e+05,1e+06)) + scale_colour_viridis() + theme(text = element_text(size = 15))
dev.off()

print("Plotting flowcell channels summary histograms")
png(filename = file.path(output.dir, "channel_summary.png"), width = 2400, height = 960) 
c = channel.summary(subset(d, Q_cutoff=="All reads"))
c10 = channel.summary(subset(d, Q_cutoff=="Reads with mean Q score > 10"))
c$Q_cutoff = "All reads"
c10$Q_cutoff = "Reads with mean Q score > 10"
cc = rbind(c, c10)
ggplot(cc, aes(x = value)) + geom_histogram(bins = 30) + facet_grid(Q_cutoff~variable, scales="free") + theme(text = element_text(size = 20))
dev.off()

