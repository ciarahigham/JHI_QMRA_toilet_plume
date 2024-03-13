##############################################################################
### Plots for the infection risk for counter A and B                       ###
### Plots infection risk for SARS-CoV-2 and norovirus                      ###
### i.e. generates 4 plots with violins for male/female                    ###
### Need to run dose_response_monte_carlo.R first for all combinations     ###
###                                                                        ###
### Ciara A. Higham March 2023 - JHI paper                                ###
##############################################################################

###############################
###        Pre-lims         ###
###############################

#install.packages("ggh4x")
#install.packages("ggplot2")

library(ggh4x)
library(ggplot2)

##############################################################################################################
##############################################################################################################

# The following blocks of code extend ggplot2 GeomViolin and make adjustments to add quantiles to the plot
# see stackoverflow (https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2) (@jan-glx, @Axeman)
# and stackoverflow (https://stackoverflow.com/questions/47651868/split-violin-plot-with-ggplot2-with-quantiles) (@Axeman)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                           data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

##############################################################################################################
##############################################################################################################

# comment out counter_a/counter_b, SARS/Norovirus

counter_loc = "counter_a"
#counter_loc = "counter_b"

pathogen = "SARS-CoV-2"
pathogen = "Norovirus"

# comparison times
t1 <- "0"
t2 <- "60"
t3 <- "240"

current_wd <- dirname(rstudioapi::getActiveDocumentContext()$path) # find current wd

# wd for dose response data from dose_response_monte_carlo.R
# and wd for location to save file
if (pathogen == "SARS-CoV-2"){
  data_wd <- paste0(current_wd,"/dose_response_output/sars_cov_2/")
  save_wd <- paste0(current_wd,"/plots/infection_risk/sars_cov_2/")
  file_name <- "sars_cov_2"
  } else if (pathogen == "Norovirus"){
    data_wd <- paste0(current_wd,"/dose_response_output/norovirus/")
    save_wd <- paste0(current_wd,"/plots/infection_risk/norovirus/")
    file_name <- "norovirus"
  }

setwd(data_wd) # set to relavant dose-response data wd

# read in the dose-response data
if (counter_loc == "counter_a"){
  male_S1 <- read.csv("male_s1_counter_a.csv",header = TRUE)
  male_S2 <- read.csv("male_s2_counter_a.csv",header = TRUE)
  female_S1 <- read.csv("female_s1_counter_a.csv",header = TRUE)
  female_S2 <- read.csv("female_s2_counter_a.csv",header = TRUE)
  file_name <- paste0(file_name, "_counter_a.png")
  
} else if (counter_loc == "counter_b"){
  male_S1 <- read.csv("male_s1_counter_b.csv",header = TRUE)
  male_S2 <- read.csv("male_s2_counter_b.csv",header = TRUE)
  female_S1 <- read.csv("female_s1_counter_b.csv",header = TRUE)
  female_S2 <- read.csv("female_s2_counter_b.csv",header = TRUE)
  file_name <- paste0(file_name, "_counter_b.png")
}

# list of the new strip labels
strip_labels_new <- c("1.5 ACH" , "3 ACH" , "6 ACH" , paste(t1, "s"), paste(t2, "s"), paste(t3, "s"), "S2", "S1")

# list of the current labels
strip_labels_current <- c("1.5", "3", "6", t1, t2, t3, "S2", "S1")

# set the strip_labels to the new strip labels
strip_labels <- setNames(strip_labels_new, strip_labels_current)

# combine all 4 sets of data
combined <- rbind(female_S1,female_S2, male_S1, male_S2)
combined$Response <- combined$Response * 100
# filter by the relevant t_enter
filtered_combined <- subset(combined, t_enter %in% c(as.numeric(t1), as.numeric(t2), as.numeric(t3)))

# set the colour for female
filtered_combined$fill_color <- ifelse(filtered_combined$Gender == "Female", "#009ADE", "#FF1F5B")

# convert the risk to a percentage


# create the plot
p <- ggplot() +
  
  # create the split violin plot
  geom_split_violin(data = filtered_combined, # select the data to plot
                    aes(x = Arrangement, y = Response, fill = fill_color), # x axis is by arrangement
                    draw_quantiles = c(0.25, 0.5, 0.75), # 25th, 50th and 75th quantile
                    alpha = 0.8, size = 0.15) +
  
  # nest the facets it is split by ACH, then split by t_enter, then split by S1/S2
  # and add the labels 
  facet_nested(. ~ ACH + t_enter + Arrangement, scales = "free", labeller = as_labeller(strip_labels)) +
  
  # plot on a log scale 
  scale_y_log10(limits = c(3.57e-10, 100), breaks = c(1e-10, 1e-8,1e-6, 1e-4,1e-2,1,100), labels=c(expression(10^{-10}) ,expression(10^{-8}), expression(10^{-6}), expression(10^{-4}), expression(10^{-2}),  "1", "100"))+
  theme_bw() +
  
  # colours for the violin
  scale_fill_manual(values = c( "#009ADE","#FF1F5B"), labels = c("Female", "Male"), name = NULL, guide = "legend") +
  scale_x_discrete(breaks = NULL) + # remove x breaks
  
  # aesthetics for the plot text size etc
  theme(
   text = element_text(family = "Arial"), # arial font
   axis.text = element_text(size = 8),
   strip.text = element_text(size = 6, margin = margin(2.65, 0.5, 2.65, 0.5)),
   strip.background = element_rect(color = "black", size = 0.2),
   legend.position = "bottom",
   legend.text = element_text(size=6),
   panel.border = element_blank(),
   axis.title.x = element_blank(),
   axis.text.x = element_blank(),
   legend.margin=margin(0.01,0.1,0.01,0.1),
   legend.key.size = unit(0.3, 'cm'),
 ) +
  
  # y label
  ylab("Risk (%)")

# print plot
print(p)

setwd(save_wd) # set save wd

# save plot
ggsave(file_name, plot = p, width = 18, height = 7, units = "cm")
