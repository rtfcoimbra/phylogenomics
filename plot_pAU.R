# Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
# data: a data frame
# measurevar: the name of a column that contains the variable to be summariezed
# groupvars: a vector containing names of columns that contain grouping variables
# na.rm: a boolean that indicates whether to ignore NA's
# conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = 0.95, .drop = TRUE) {
  library(plyr)
  # version of length() which can handle NA's: if na.rm == T, don't count them
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      length(x)
    }
  }
  # for each group's data frame, return a vector with N, mean, and sd
  datac <- ddply(data, groupvars, .drop = .drop,
                 .fun = function(xx, col) {
                   c(N = length2(xx[[col]], na.rm = na.rm),
                     mean = mean(xx[[col]], na.rm = na.rm),
                     median = median(xx[[col]], na.rm = na.rm),
                     sd = sd(xx[[col]], na.rm = na.rm))
                   },
                 measurevar)
  # rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  # calculate standard error of the mean
  datac$se <- datac$sd / sqrt(datac$N)
  # confidence interval multiplier for standard error
  # calculate t-statistic for confidence interval:
  # e.g., if conf.interval is 0.95, use 0.975 (above/below), and use df = N - 1
  ciMult <- qt(conf.interval / 2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

# load libraries
library(ggplot2)
library(ggtree)

# set up working directory
setwd("~/Sync/rcoimbra_phd/results")

# read file with tested tree topologies
trees <- read.tree("./au_test/topologies.tree")
# plot topologies
ggtree(trees, ladderize = TRUE, right = TRUE) +
  facet_wrap(~.id, ncol = 3) +
  geom_tiplab(size = 2.75) +
  ggplot2::xlim(0, 6.25) #+ theme(strip.background = element_blank())
# save plot in '.pdf' format
ggsave("./au_test/topologies.pdf", width = 170, height = 257, units = "mm")

# read input
df <- read.table("./au_test/combined.au", sep = "\t", header = TRUE)
# change column names
colnames(df) <- c("frag_size", "id", "pAU")
# remove rows with missing data from data frame
df <- df[complete.cases(df), ]
# obtain summary statistics with custom function
df.s <- summarySE(df, measurevar = "pAU", groupvars = c("frag_size", "id"))

# show when significant rejection is reached
head(subset(df.s, id == "Top1" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top2" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top3" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top4" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top5" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top6" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top7" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top8" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top9" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top10" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top11" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top12" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top13" & pAU - ci > 0.95))
head(subset(df.s, id == "Top14" & pAU + ci <= 0.05))
head(subset(df.s, id == "Top15" & pAU + ci <= 0.05))

# plot pAU over fragment sizes
ggplot(df.s, aes(color = id, x = frag_size, y = pAU)) +
  geom_point(size = 1.1) +
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin = pAU - ci, ymax = pAU + ci), width = 0.05) +
  geom_hline(aes(yintercept = 0.95), color = "green", lty = 2) +
  geom_hline(aes(yintercept = 0.05), color = "red", lty = 2) +
  ylab("AU values") +
  xlab("Fragment size") +
  scale_color_viridis_d() +
  scale_x_discrete(limit = c(seq(100000, 600000, 100000)),
                   labels = c("100", "200", "300", "400", "500", "600")) +
  theme_classic() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 7),
        legend.title = element_blank())
# save plot in '.pdf' format
ggsave("./au_test/au_test.pdf", width = 170, height = 85, units = "mm")
