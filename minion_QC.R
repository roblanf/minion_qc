args<-commandArgs(TRUE)

input.file = args[1]
output.dir = args[2]

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
    return(length(which(d$sequence_length_template>=len)))
}

bases.gt <- function(d, len){
  reads = subset(d, sequence_length_template >= len)
  return(sum(as.numeric(reads$sequence_length_template)))
}

summary.stats <- function(d){
    
    total.bases = sum(as.numeric(d$sequence_length_template))
    N50.length = d$sequence_length_template[min(which(d$cumulative.bases > (total.bases/2)))]
    mean.length = mean(as.numeric(d$sequence_length_template))
    median.length = median(as.numeric(d$sequence_length_template))
    max.length = max(as.numeric(d$sequence_length_template))
    mean.q = mean(d$mean_qscore_template)
    median.q = median(d$mean_qscore_template)
    
    reads = c(reads.gt(d, 20000), 
              reads.gt(d, 50000),
              reads.gt(d, 100000),
              reads.gt(d, 200000),
              reads.gt(d, 500000),
              reads.gt(d, 1000000))
    names(reads) = c(">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m")

    bases = c(bases.gt(d, 20000), 
              bases.gt(d, 50000),
              bases.gt(d, 100000),
              bases.gt(d, 200000),
              bases.gt(d, 500000),
              bases.gt(d, 1000000))
    names(bases) = c(">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m")
    
                
    return(list('total.bases' = total.bases, 
                'N50.length' = N50.length, 
                'mean.length' = mean.length, 
                'median.length' = median.length,
                'max.length' = max.length,
                'mean.q' = mean.q,
                'median.q' = median.q,
                'reads' = reads,
                'bases' = bases
                ))
}

channel.summary <- function(d){
    a = ddply(d, .(channel), summarize, total.bases = sum(sequence_length_template), total.reads = sum(which(sequence_length_template>=0)), mean.read.length = mean(sequence_length_template), median.read.length = median(sequence_length_template))
    b = melt(a, id.vars = c("channel"))
    return(b)    
}

fnlist <- function(x, fil){ z <- deparse(substitute(x))
cat(z, "\n", file=fil, append=TRUE)
nams=names(x) 
for (i in seq_along(x) ){ 
  cat("", "\t", paste(names(x[[i]]), "\t"), "\n", file=fil, append=TRUE)
  cat(nams[i], "\t",  x[[i]], "\n", file=fil, append=TRUE) }
}



# supress warnings
options(warn=-1)


print("Creating output directory")
dir.create(output.dir)
out.txt = file.path(output.dir, "summary.txt")

write(paste("Summary stats from input file", input.file, "\n\n"), file = out.txt)


# write summaries
print("Loading and summarising input file")
d = load_summary(input.file)
all.reads.summary = summary.stats(d)
fnlist(all.reads.summary, out.txt)

d10 = load_summary(input.file, min.q = 10)
write("\n\n", file = out.txt, append=TRUE)
q10.reads.summary = summary.stats(d10)
fnlist(q10.reads.summary, out.txt)

d$Q_cutoff = "All reads"
d10$Q_cutoff = "Reads with mean Q score > 10"
dd = rbind(d, d10)


# make plots
print("Plotting length histogram")
png(filename = file.path(output.dir, "length_histogram.png"), width = 960, height = 960)
ggplot(dd, aes(x = sequence_length_template)) + geom_histogram(bins = 500) + scale_x_log10() + facet_wrap(~Q_cutoff, ncol = 1)
dev.off()

print("Plotting mean Q score histogram")
png(filename = file.path(output.dir, "q_histogram.png"), width = 960, height = 960)
ggplot(dd, aes(x = mean_qscore_template)) + geom_histogram(bins = 500) + facet_wrap(~Q_cutoff, ncol = 1)
dev.off()

print("Plotting events per base histogram")
dd$events_per_base = dd$num_events_template/dd$sequence_length_template
png(filename = file.path(output.dir, "epb_histogram.png"), width = 960, height = 960)
ggplot(dd, aes(x = events_per_base)) + geom_histogram(bins = 500) + scale_x_log10() + facet_wrap(~Q_cutoff, ncol = 1)
dev.off()

print("Plotting read length vs. q score scatterplot")
d$events_per_base = d$num_events_template/d$sequence_length_template
d$events_per_base[which(d$events_per_base>10)] = 10
png(filename = file.path(output.dir, "length_vs_q.png"), width = 960, height = 960)
ggplot(d, aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + geom_point(alpha=0.05, size = 0.4) + scale_x_log10(breaks=c(1e+01, 1e+02, 1e+03,1e+04,1e+05,1e+06)) + scale_colour_viridis()
dev.off()

print("Plotting flowcell channels events per base plot")
d$events_per_base = d$num_events_template/d$sequence_length_template
png(filename = file.path(output.dir, "flowcell_channels_epb.png"), width = 2400, height = 2400)
ggplot(d, aes(x=start_time/3600, y=events_per_base, colour = mean_qscore_template)) + 
  geom_point(size=1, alpha=0.5) + 
  scale_colour_viridis() + 
  scale_y_log10(limits=c(1e+0, 1e+01)) + 
  facet_grid(row~col) +
  theme(panel.spacing = unit(0.5, "lines")) +
  xlab("Hours into run") +
  theme(legend.position="none")
dev.off()

print("Plotting flowcell channels summary histograms")
png(filename = file.path(output.dir, "channel_summary.png"), width = 2400, height = 960)
c = channel.summary(d)
c10 = channel.summary(d10)
c$Q_cutoff = "All reads"
c10$Q_cutoff = "Reads with mean Q score > 10"
cc = rbind(c, c10)
ggplot(cc, aes(x = value)) + geom_histogram(bins = 30) + facet_grid(Q_cutoff~variable, scales="free")
dev.off()

print("Plotting flowcell yield summary")
png(filename = file.path(output.dir, "yield_summary.png"), width = 960, height = 960)
ggplot(dd, aes(x=sequence_length_template, y=cumulative.bases, colour = Q_cutoff)) + 
  geom_line(size = 1) + 
  scale_x_continuous(breaks =c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000), limits = c(0, 100000)) +
  xlab("Minimum read length") +
  ylab("Total yield in bases")
dev.off()
