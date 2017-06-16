# look at albacore's sequencing summary file
library(ggplot2)
library(viridis)
library(plyr)
library(reshape2)

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


load_summary <- function(filepath, flowcell="NA", min.q=0){
    # load a sequencing summary and add some info
    d = read.delim(filepath)
    d$sequence_length_template = as.numeric(d$sequence_length_template)
    d = subset(d, mean_qscore_template >= min.q)
    d = merge(d, map, by="channel")
    d$flowcell = flowcell
    d = d[with(d, order(-sequence_length_template)), ] # sort by read length
    d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))
    return(d)
}

reads.gt <- function(d, len){
    return(length(which(d$sequence_length_template>len)))
}


summary.stats <- function(d){
    
    total.bases = sum(as.numeric(d$sequence_length_template))
    N50.length = d$sequence_length_template[min(which(d$cumulative.bases > (total.bases/2)))]
    mean.length = mean(as.numeric(d$sequence_length_template))
    median.length = median(as.numeric(d$sequence_length_template))
    mean.q = mean(d$mean_qscore_template)
    median.q = median(d$mean_qscore_template)
    
    reads = c(reads.gt(d, 20000), 
              reads.gt(d, 50000),
              reads.gt(d, 100000),
              reads.gt(d, 200000),
              reads.gt(d, 500000),
              reads.gt(d, 1000000))
    names(reads) = c(">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m")
            
    return(list('total.bases' = total.bases, 
                'N50.length' = N50.length, 
                'mean.length' = mean.length, 
                'median.length' = median.length,
                'mean.q' = mean.q,
                'median.q' = median.q,
                'reads' = reads
                ))
}

channel.summary <- function(d){
    a = ddply(d, .(channel), summarize, total.bases = sum(sequence_length_template), total.reads = sum(which(sequence_length_template>=0)), mean.read.length = mean(sequence_length_template), median.read.length = median(sequence_length_template))
    b = melt(a, id.vars = c("channel"))
    p = ggplot(b, aes(x = value)) + geom_histogram() + facet_wrap(~variable, scales="free", ncol = 1)
    return(list("channels" = b, "plot" = p))    
}

d = load_summary("Desktop/summaries/sequencing_summary_A2.txt")
summary.stats(d)

ggplot(d, aes(x = sequence_length_template)) + geom_histogram() + scale_x_log10(breaks=c(1e+01, 1e+02, 1e+03,1e+04,1e+05,1e+06))
ggplot(d, aes(x = mean_qscore_template)) + geom_histogram()
ggplot(d, aes(x = sequence_length_template, y = mean_qscore_template)) + geom_point(alpha=0.01, size = 0.1) + scale_x_log10(breaks=c(1e+01, 1e+02, 1e+03,1e+04,1e+05,1e+06))
d$events_per_base = d$num_events_template/d$sequence_length_template
ggplot(d, aes(x = events_per_base)) + geom_histogram() + scale_x_log10()
d$events_per_base[which(d$events_per_base>10)] = 10
ggplot(d, aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + geom_point(alpha=0.01, size = 0.1) + scale_x_log10(breaks=c(1e+01, 1e+02, 1e+03,1e+04,1e+05,1e+06)) + scale_colour_viridis()
d$events_per_base = d$num_events_template/d$sequence_length_template
a = subset(d, channel == 213)
ggplot(a, aes(x = start_time/3600, y = events_per_base, colour = mean_qscore_template)) + geom_point(size = 1, alpha = 0.5) + xlab("Hours into run") + scale_y_log10() + scale_colour_viridis() + xlim(0, 35)
a = subset(d, channel == 408)
ggplot(a, aes(x = start_time/3600, y = events_per_base, colour = mean_qscore_template)) + geom_point(size = 1, alpha = 0.5) + xlab("Hours into run") + scale_y_log10() + scale_colour_viridis() + xlim(0, 35)

ggplot(d, aes(x=start_time/3600, y=events_per_base, colour = mean_qscore_template)) + 
  geom_point(size=0.1, alpha=0.2) + 
  scale_colour_viridis() + 
  scale_y_log10(limits=c(1e+0, 1e+01)) + 
  facet_grid(row~col) +
  theme(panel.spacing = unit(0.1, "lines")) +
  xlab("Hours into run") +
  theme(legend.position="none")

channel.summary(d)

#### 

d10 = load_summary("Desktop/summaries/sequencing_summaryG1.txt", min.q = 10)
summary.stats(d10)


# cumulative bases by read length
ggplot(d, aes(x=sequence_length_template, y=cumulative.bases)) + geom_line() + scale_x_log10() + scale_y_log10()
ggplot(d, aes(x=sequence_length_template, y=cumulative.bases)) + geom_line() + scale_x_continuous(breaks =c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000), limits = c(0, 100000))


# e.g. to compare two flowcells
f1 = load_summary("Desktop/summaries/sequencing_summary_A2.txt", flowcell = "A2", min.q = 0)
f2 = load_summary("Desktop/summaries/sequencing_summaryG1.txt", flowcell = "G1", min.q = 0)
f = rbind(f1, f2)

f3 = load_summary("Desktop/summaries/sequencing_summary_A2.txt", flowcell = "A2", min.q = 10)
f4 = load_summary("Desktop/summaries/sequencing_summaryG1.txt", flowcell = "G1", min.q = 10)
fq = rbind(f3, f4)

f$min.q = "minimum mean qscore of reads: 0"
fq$min.q = "minimum mean qscore of reads: 10"
fd = rbind(f, fq)

# yield from reads at least as long as x
ggplot(fd, aes(x=sequence_length_template, y=cumulative.bases, colour = flowcell)) + 
    geom_line() + 
    scale_x_continuous(breaks =c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000), limits = c(0, 100000)) +
    facet_wrap(~min.q, ncol=1) +
    xlab("Minimum read length") +
    ylab("Total yield in bases")
    
ggplot(f, aes(x=mean_qscore_template, colour = flowcell)) + geom_density(aes(y = ..count..))
ggplot(f, aes(x=sequence_length_template, colour = flowcell)) + geom_density(aes(y = ..count..)) + scale_x_log10(breaks = c(100, 1000, 10000, 20000, 50000, 100000)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(f, aes(x=sequence_length_template, colour = flowcell)) + geom_density(aes(y = ..count..)) + xlim(0, 100000)
ggplot(f, aes(x=sequence_length_template, colour = flowcell)) + geom_density(aes(y = ..count..)) + xlim(10000, 100000)

ggplot(fd, aes(x=mean_qscore_template, colour = flowcell)) + geom_density(aes(y=..count..))
ggplot(fd, aes(x=sequence_length_template, colour = flowcell)) + geom_density(aes(y = ..count..)) + scale_x_log10(breaks = c(100, 1000, 10000, 20000, 50000, 100000))
ggplot(fd, aes(x=sequence_length_template, colour = flowcell)) + geom_density(aes(y = ..count..)) + xlim(0, 100000)

ggplot(fd, aes(x=sequence_length_template, y=mean_qscore_template)) + 
    stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) + 
    scale_fill_viridis() + facet_wrap(~flowcell, ncol=1) + xlim(0, 50000)

s1 = channel.summary(f1)
s2 = channel.summary(f2)
s1$channels$flowcell = "A2"
s2$channels$flowcell = "G1"
s = rbind(s1$channels, s2$channels)
ggplot(s, aes(x = value)) + geom_histogram(aes(fill = flowcell)) + facet_grid(flowcell~variable, scales="free") + ggtitle("Stats per channel for 2 flowcells")
