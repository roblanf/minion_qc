# look at albacore's sequencing summary file
library(ggplot2)
library(viridis)


load_summary <- function(filepath, flowcell="NA", min.q=0){
    # load a sequencing summary and add some info
    d = read.delim(filepath)
    d = subset(d, mean_qscore_template >= min.q)
    map = data.frame(channel=1:512, row=rep(1:16, 32), col=rep(1:32, each=16))   
    d = merge(d, map, by="channel")
    d$flowcell = flowcell
    d = d[with(d, order(-sequence_length_template)), ] # sort by read length
    d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))
    return(d)
}


d = load_summary("Desktop/summaries/sequencing_summaryG1.txt")

ggplot(d, aes(x=start_time/3600, y=sequence_length_template, colour = mean_qscore_template)) + geom_point(size=0.2, alpha=0.5) + scale_y_log10(labels = scales::comma, breaks = c(0, 50, 100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 750000, 1000000)) + scale_colour_viridis()

ggplot(d, aes(x=start_time/3600, y=num_events, colour = mean_qscore_template)) + geom_point(size=0.2, alpha=0.5) + scale_y_log10(labels = scales::comma) + scale_colour_viridis()

ggplot(d, aes(x=start_time/3600, y=num_events/sequence_length_template, colour = mean_qscore_template)) + 
  geom_point(size=0.2, alpha=0.5) + 
  scale_y_log10(labels = scales::comma) + 
  scale_colour_viridis()

ggplot(d, aes(x=sequence_length_template, y=num_events, colour = mean_qscore_template)) + 
  geom_point(size=0.2, alpha=0.2) + 
  scale_x_log10(labels = scales::comma, breaks = c(10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_y_log10(labels = scales::comma, breaks = c(10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_colour_viridis() + 
  geom_abline(slope = 1, gradient = 1)

ggplot(d, aes(x=sequence_length_template, y=num_events/sequence_length_template, colour = mean_qscore_template)) + 
  geom_point(size=0.2, alpha=0.2) + 
  scale_x_log10(labels = scales::comma, breaks = c(10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_y_log10(labels = scales::comma, breaks = c(10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_colour_viridis() + 
  geom_hline(yintercept = 1, linetype='dotted', colour='red')

ggplot(d, aes(x=sequence_length_template, y=num_events/sequence_length_template, colour = mean_qscore_template)) + 
  geom_point(size=0.2, alpha=0.2) + 
  scale_x_log10(labels = scales::comma, breaks = c(1, 2, 5, 10, 100, 1000, 10000, 100000, 1000000), limits=c(0, 10)) + 
  scale_y_log10(labels = scales::comma, breaks = c(1, 2, 5, 10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_colour_viridis() + 
  geom_hline(yintercept = 1, linetype='dotted', colour='red') +
  facet_grid(row~col) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank(), panel.spacing = unit(0.01, "lines"))

ggplot(d, aes(x=start_time/3600, y=log(num_events/sequence_length_template), colour = mean_qscore_template)) + 
  geom_point(size=0.2, alpha=0.2) + 
  scale_colour_viridis() + 
  ylim(0,4) + 
  facet_wrap(~channel, ncol = 16) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank(), panel.spacing = unit(0.1, "lines"))

ggplot(d, aes(x=sequence_length_template)) + geom_histogram(aes(y=cumsum(..count..)), bins = 100) + scale_x_log10() 

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

