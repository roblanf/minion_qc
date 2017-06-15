# look at albacore's sequencing summary file
library(ggplot2)
library(viridis)

d = read.delim("Desktop/sequencing_summary.txt")
d = d[with(d, order(sequence_length_template)), ] # sort by read length
d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))

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
  scale_x_log10(labels = scales::comma, breaks = c(1, 2, 5, 10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_y_log10(labels = scales::comma, breaks = c(1, 2, 5, 10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_colour_viridis() + 
  geom_hline(yintercept = 1, linetype='dotted', colour='red') +
  facet_wrap(~channel) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank(), panel.spacing = unit(0.01, "lines"))

ggplot(d, aes(x=start_time/3600, y=log(num_events/sequence_length_template), colour = mean_qscore_template)) + 
  geom_point(size=0.2, alpha=0.2) + 
  scale_colour_viridis() + 
  ylim(0,4) + 
  facet_wrap(~channel, ncol = 16) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank(), panel.spacing = unit(0.1, "lines"))

ggplot(d, aes(x=sequence_length_template)) + geom_histogram(aes(y=cumsum(..count..)), bins = 100) + scale_x_log10() 

# cumulative bases by read length
ggplot(d, aes(x=sequence_length_template, y=cumulative.bases)) + geom_line() + scale_x_log10()
ggplot(d, aes(x=sequence_length_template, y=cumulative.bases)) + geom_line() + xlim(0, 100000)

