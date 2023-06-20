# Packages
pacman::p_load(here, 
               DESeq2, # BiocManager::install("DESeq2")
               tximport, # BiocManager::install("tximport")
               rhdf5, # BiocManager::install("DESeq2")
               microseq, 
               rtracklayer, # BiocManager::install("rtracklayer")
               MASS, 
               tidyverse, 
               patchwork, 
               ggrepel,
               ggpubr)

#TODO: 
  #FIGS2
  #FIGS3
  #FIGS4


#### Global variables ####
lfc_thresh = 2
p_thresh = 0.01
p_adj_meth = "fdr"
custom_palette = RColorBrewer::brewer.pal(8, "Set2")
brew_blue = RColorBrewer::brewer.pal(3, "Blues")[2]

custom_red = "#c7546e"
custom_grey = "#B3B3B3"
cust_mag = "#E78AC3"
cust_green = "#A6D854"
text_size = 10
  
#### functions ####
ratio_plot_fill <- function(df, 
                            mapping, 
                            text_pos,
                            strain = NULL) {
  if (!is.null(strain)) {
    reference <- reference %>% 
      dplyr::filter(str_detect(var, {{strain}}))
  }
  
  ggplot(df, mapping) +
    geom_hline(data = dplyr::filter(reference, interaction == "bsps"),
               aes(yintercept = mean),
               linetype = 2) +
    geom_rect(data = dplyr::filter(reference, interaction == "bsps"),
              aes(ymax = max,
                  ymin = min),
              xmin = -Inf, 
              xmax = Inf,
              alpha = .1,
              fill = "black",
              inherit.aes = F) +
    geom_hline(data = dplyr::filter(reference, interaction == "dhbps"),
               aes(yintercept = mean),
               linetype = 2,
               col = "#A90F33") +
    geom_rect(data = dplyr::filter(reference, interaction == "dhbps"),
              aes(ymax = max,
                  ymin = min),
              xmin = -Inf, 
              xmax = Inf,
              alpha = .1,
              fill = "#A90F33",
              inherit.aes = F) +
    geom_text(data = dplyr::filter(reference, interaction == "bsps"),
              aes(y = ifelse(max < 1, 
                             max - 0.3, 
                             max - 0.04)),
              x = text_pos + 0.5,
              label = "Bs vs Ps",
              hjust = 1,
              size = 5,
              fontface = "bold",
              inherit.aes = F) +
    geom_text(data = dplyr::filter(reference, interaction == "dhbps"),
              aes(y = max - 0.04),
              x = text_pos + 0.5,
              label = "dhbA vs Ps",
              hjust = 1,
              size = 5,
              fontface = "bold",
              col = "#A90F33",
              inherit.aes = F) +
    geom_boxplot(outlier.shape = NA) +
    ggforce::geom_sina() +
    scale_y_continuous(limits = c(0.1, 1.7),
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)) +
    theme_bw() +
    labs(y = "Area coculture / area monoculture") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = text_size),
      text = element_text(size = text_size)
    )
}

ratio_plot_gradient <- function(df, 
                                mapping, 
                                text_pos,
                                strain = NULL) {
  if (!is.null(strain)) {
    reference <- reference %>% 
      dplyr::filter(str_detect(var, {{strain}}))
  }
  
  ggplot(df, mapping) +
    geom_hline(data = dplyr::filter(reference, interaction == "bsps"),
               aes(yintercept = mean),
               linetype = 2) +
    geom_rect(data = dplyr::filter(reference, interaction == "bsps"),
              aes(ymax = max,
                  ymin = min),
              xmin = -Inf, 
              xmax = Inf,
              alpha = .1,
              fill = "black",
              inherit.aes = F) +
    geom_hline(data = dplyr::filter(reference, interaction == "dhbps"),
               aes(yintercept = mean),
               linetype = 2,
               col = "#A90F33") +
    geom_rect(data = dplyr::filter(reference, interaction == "dhbps"),
              aes(ymax = max,
                  ymin = min),
              xmin = -Inf, 
              xmax = Inf,
              alpha = .1,
              fill = "#A90F33",
              inherit.aes = F) +
    geom_text(data = dplyr::filter(reference, interaction == "bsps"),
              aes(y = ifelse(max < 1, 
                             max - 0.3, 
                             max - 0.04)),
              x = text_pos + 0.5,
              label = "Bs vs Ps",
              hjust = 1,
              size = 4,
              fontface = "bold",
              inherit.aes = F) +
    geom_text(data = dplyr::filter(reference, interaction == "dhbps"),
              aes(y = max - 0.04),
              x = text_pos + 0.5,
              label = "dhbA vs Ps",
              hjust = 1,
              size = 4,
              fontface = "bold",
              col = "#A90F33",
              inherit.aes = F) +
    geom_boxplot(outlier.shape = NA,
                 fatten = 1) +
    ggforce::geom_sina(alpha = 0.5,
                       size = 1) +
    scale_y_continuous(limits = c(0.1, 1.7),
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)) +
    theme_bw() +
    labs(y = "Area coculture / area monoculture") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = text_size),
      text = element_text(size = text_size)
    )
}

parse_plate <- function(data, 
                        plate_layout, 
                        time_series = F, 
                        wide_data = T, 
                        wide_layout = T) {
  ## Check if layout has been provided in wide format. If yes, convert to long format
  if (wide_layout == T) {
    layout_names_idx <- grep(
      pattern = "[A-Za-z]{2}", # Looking for cells with at least two letters
      plate_layout[, 1]
    )
    wells <- as.character(plate_layout[grep(
      pattern = "^[A-Z]$",
      plate_layout[, 1]
    ), 1])
    
    ## Create data frame to store information
    ## A bit clunky, though it is able to handle all types of titer plates
    layout_long_test <-
      data.frame(well = c(paste(LETTERS[seq(
        from = 1,
        to = max(which(LETTERS[1:26] %in% wells))
      )],
      rep(1:(ncol(plate_layout) - 1),
          each = max(which(LETTERS[1:26] %in% wells))
      ),
      sep = ""
      )))
    
    ## The layout file block should start at the row right after the variable name
    for (i in seq_len(length(layout_names_idx))) {
      block_name <- plate_layout[layout_names_idx[i], 1]
      
      block_start <- layout_names_idx[i] + 1
      block_end <- layout_names_idx[i] + length(unique(wells))
      
      new_block <- plate_layout[(block_start):(block_end), 2:ncol(plate_layout)]
      new_block <- reshape2::melt(new_block, measure.vars = names(new_block))
      
      layout_long_test[, block_name] <- new_block[, 2]
    }
    plate_layout <- layout_long_test
  }
  
  ## If data is end point, find measurement blocks from raw file information
  ## This only works for the tecan plate reader, so far. I don't know what raw files from the
  ## Synergy looks like yet, so I can't grep for anything.
  if (time_series == F) {
    start_time_idx <- which(data[, 1] == "Start Time:")
    end_idx <- which(data[, 1] == "End Time:")
    names_idx <- grep(pattern = "Label:.*", data[, 1])
    
    all_data <- c()
    for (i in seq_len(length(start_time_idx))) {
      block_name <- data[names_idx[i], 1]
      block_name <- gsub(pattern = "Label: ", replacement = "", block_name)
      
      block_start <- start_time_idx[i] + 4
      block_end_idx <- end_idx[i] - 5
      
      ## If data was provided in wide format, it must be reformatted into a long format
      if (wide_data == T) {
        new_block <- data[(block_start):(block_end_idx), 1:7]
        new_block <- reshape2::melt(new_block, measure.vars = names(new_block[2:7]))
        new_block[, 2] <- rep(x = c(1:6), each = 4)
        new_block <- transform(new_block, well = paste(V1, variable, sep = ""))
        new_block <- cbind(new_block[4], new_block[3])
        
        ## Else if wide_data == F
      } else {
        new_block <- data[(block_start):(block_end_idx), 1:2]
        names(new_block)[1] <- "well"
        names(new_block)[2] <- "value"
      }
      
      ## Join parsed data to layout data to add layout variables
      new_block$value <- as.numeric(new_block$value)
      joined_block <- dplyr::full_join(plate_layout, new_block, by = "well")
      joined_block$measure <- block_name
      
      all_data <- rbind(all_data, joined_block)
    }
  }
  ## Else if time_series == T, adjust measurement blocks accordingly
  ## This only works for the Synergy HTX.
  ## I don't know what time series data from the tecan looks like, yet.
  else {
    names_idx <- which(grepl(pattern = "[1-9]:.*", data[, 1])) # Looking for identifiers like '1:...'
    all_data <- c()
    for (i in seq_len(length(names_idx))) {
      block_start <- names_idx[i] + 2
      block_end <- names_idx[i] + min(which(data[block_start:nrow(data), 3] %in% c("", NA)))
      block_name <- sub(pattern = "Read [1-9]:", "", data[names_idx[i], 1])
      
      new_block <- data[(block_start):(block_end), 2:ncol(data)]
      names(new_block) <- new_block[1, ]
      new_block <- new_block[-1, -2]
      new_block <- reshape2::melt(new_block, measure.vars = names(new_block[-1]), variable.name = "well")
      
      ## Join parsed data to layout data to add layout variables
      joined_block <- dplyr::full_join(plate_layout, new_block, by = "well")
      joined_block$measure <- block_name
      all_data <- rbind(all_data, joined_block)
    }
  }
  
  if (any(is.na(all_data$value))) warning("You have missing values in your data. Are you sure you provided the right format?")
  return(all_data)
} # End function


# LCMS parsing function. Takes a raw masshunter file read in as a data.frame
# e.g. df <- read.csv(x, header = F)
# Also takes chr_types (a character vector with the chromatograms in the data; ESI or DAD)
parse_chromatograms <- function(chromatogram_csv, chr_types) {
  block_start <- which(str_detect(
    chromatogram_csv$V1,
    str_c(chr_types,
          collapse = "|"
    )
  ))
  block_end <- c(
    block_start[-1] - 1,
    nrow(chromatogram_csv)
  )
  
  blocks <- map2(block_start, block_end, function(start, end) {
    block <- chromatogram_csv[start:end, ]
    names(block) <- block[2, ]
    block <- block[-2, ] %>%
      mutate(info = .[1, 1]) %>%
      as_tibble() %>%
      dplyr::select(-1)
    block <- block[-1, ] %>%
      rename(
        time_m = "X(Minutes)",
        value = 2
      ) %>%
      mutate(
        time_m = as.numeric(time_m),
        value = as.numeric(value)
      )
  })
  
  blocks_parsed <- blocks %>%
    map(function(chromatogram) {
      map(chr_types, function(type) {
        if (type == "ESI") {
          if (all(str_detect(chromatogram$info, "ESI"))) {
            chromatogram %>% 
              separate(info,
                       into = c("mode", "chromatogram", "ms_level", NA, "file"),
                       sep = " ")
          }
        } else if (type == "DAD") {
          if (all(str_detect(chromatogram$info, "DAD"))) {
            chromatogram %>%
              separate(info,
                       into = c("mode", NA, "range", "file"),
                       sep = " ")
          }
        }
      }) %>% 
        set_names(chr_types) %>% 
        compact() 
    }) %>% 
    list_flatten(name_repair = "unique")
}


#### Figure 1: 3610-PS92 antagonism is dhbA-dependent ####
# Colony area and cfu
# Data import and manipulation
data <- map(list.files(here("data", "area", "base_interaction"), ".csv", full.names = T), read.csv) %>% 
  map2(list.files(here("data", "area", "base_interaction"), ".csv"), function(x, y) {
    x %>% 
      mutate(time = str_remove(y, "_area.csv"),
             Label = str_remove(Label, "-[:digit:]{2}"),
             Label = str_remove(Label, ".czi")) %>% 
      separate(Label, into = c("interaction", "type", "strain")) %>% 
      dplyr::select(-X) %>% 
      rename(area_um2 = Area)
  })

df <- do.call("rbind", data) %>% 
  mutate(strain = ifelse(type == "co", 
                         strain,
                         ifelse(type == "bs", "bs", "ps")),
         type = ifelse(type == "co", type, "mono"),
         area_mm2 = area_um2/10^6,
         var = str_c(strain, type, sep = "_"),
         interaction = as.factor(interaction))

# Split on type to calculate ratio
split_list <- df %>% 
  group_by(type) %>% 
  group_split()

ratio <- split_list[[2]] %>%
  group_by(interaction, strain, time) %>% 
  summarize(mono_mean = mean(area_mm2), .groups = "drop") %>% 
  right_join(ungroup(split_list[[1]]), by = c("interaction", "strain", "time")) %>% 
  mutate(l2fc = log2(area_mm2 / mono_mean),
         ratio = area_mm2/mono_mean)

# Reference for other figures
reference <- ratio %>% 
  group_by(interaction, strain, time, type, var) %>% 
  summarise(mean = mean(ratio),
            min = min(ratio),
            max = max(ratio),
            se = sd(ratio)/sqrt(length(ratio)), 
            CI = qt(0.975, df = length(ratio)-1) * sd(ratio) / sqrt(length(ratio))) %>% 
  dplyr::filter(time == "72h") %>% 
  mutate(strain = as.factor(strain))


levels(reference$strain) = c("bs" = expression(paste(italic(B.~subtilis))),
                             "ps" = expression(paste(italic(P.~marginalis)))
)

levels(ratio$interaction) = c("bsps" = "3610 + PS92",
                              "dhbps" = "dhbA + PS92"
)

  # CFU from colonies
cfu <- read_csv2(here("data", "cfu", "colony_cfu.csv")) %>% 
  mutate(cfu = count * 10^dil_factor,
         log_cfu = log10(cfu),
         interaction = as.factor(interaction))


levels(cfu$interaction) = c("bsps" = "3610 + PS92",
                            "dhbps" = "dhbA + PS92"
)

# Statistical test
cfu_test <- compare_means(count ~ type, 
                          data = cfu, 
                          method = "t.test", 
                          group.by = c("species", "interaction")) %>% 
  mutate(p.signif = ifelse(p.adj > 0.05, "ns", p.signif))


# LCMS
l <- map(list.files(here("data/lcms"), pattern = "SAL.*CSV", full.names = T), function(x) {
  df <- read.csv(x, header = F) %>% 
    parse_chromatograms(chr_types = c("ESI", "DAD"))
}) %>% set_names(list.files(here("data/lcms"), pattern = "SAL.*CSV")) # I set the names of my list elements to the file-names

md <- readxl::read_xlsx(here("data/lcms/md.xlsx"), col_names = T)

# Ionization chromatograms
esi <- map(l, function(file) {
  file[str_detect(names(file), "ESI")] %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% # Flatten to single df
  mutate(file = str_remove(file, ".d")) %>% 
  left_join(md, by = "file") %>% # Add md
  mutate(sample = fct_relevel(sample, "kb_solid", "lb_solid", "kb_liquid", "lb_liquid"))

# UV/VIS chromatograms
dad <- map(l, function(file) {
  file[str_detect(names(file), "DAD")] %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% # Flatten to single df
  mutate(file = str_remove(file, ".d")) %>% 
  left_join(md, by = "file") %>% # Add md
  mutate(sample = fct_relevel(sample, "kb_solid", "lb_solid", "kb_liquid", "lb_liquid"))

# Area under the curve (manually extracted from MassHunter)
auc <- read.table(here("data/lcms", "auc.csv"), sep = ";", dec = ".", header = T) %>%
  mutate(sample = fct_relevel(sample, "kb_solid", "lb_solid", "kb_liquid", "lb_liquid"),
         compound = fct_relevel(compound, "bb", "bb_b", "bb_c", "pulch", "bacillaene", "srf"))

# Plots
# Colony ratio
ratio_p <- ratio %>% 
  ggplot(aes(x = strain, y = ratio, fill = time)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_violin() +
  ggforce::geom_sina(alpha = 0.5) +
  stat_compare_means(inherit.aes = T,
                     method = "anova",
                     label = "p.format",
                     paired = F) +
  facet_grid(. ~ interaction, 
             labeller = label_parsed) +
  scale_fill_brewer(palette = "Blues",
                    name = "Time",
                    labels = c("24h", "48h", "72h")) +
  scale_x_discrete(labels = c("Bacillus",
                              "Pseudomonas"
                              )
  ) +
  scale_y_continuous(breaks = c(0.4, 0.8, 1.0, 1.2, 1.6)) +
  theme_bw() +
  labs(y = "Area coculture / area monoculture",
       x = "") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "bottom",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-10,-10,-5,-10)
  )

# CFU from colonies
cfu_p <- cfu %>% 
  ggplot(aes(x = species, y = cfu, fill = type)) +
  geom_violin() +
  ggforce::geom_sina(alpha = 0.5) +
  geom_pwc(tip.length = 0.01,
           method = "wilcox_test",
           label = "p.adj.signif",
           p.adjust.method = "BH",
           p.adjust.by = "group",
           hide.ns = F
  ) +
  facet_grid(. ~ interaction,
             labeller = labeller(interaction = label_parsed)) +
  scale_fill_manual(name = "Interaction type",
                    labels = c("Coculture", "Monoculture"),
                    values = c(custom_red, brew_blue)) +
  scale_x_discrete(labels = c("Bacillus",
                              "Pseudomonas"
                              )
  ) +
  scale_y_log10() +
  theme_bw() +
  labs(y = "CFU per macrocolony",
       x = "") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "bottom",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-10,-10,-5,-10)
  )

# UV/VIS LCMS
dad_p <- dad %>% 
  dplyr::filter(str_detect(sample, "solid")) %>% 
  ggplot(aes(x = time_m, 
             y = value, 
             fill = fct_reorder(rep, value, .desc = F),
             col = fct_reorder(rep, value, .desc = F))) +
  annotate(geom = "rect", 
           xmin = 4.1, 
           xmax = 4.5, 
           ymin = 0,
           ymax = Inf,
           fill = custom_red, colour = custom_red, alpha = 0.1) + 
  annotate(geom = "text",
           x = 4.65,
           y = 210,
           label = "BB") +
  annotate(geom = "rect", 
           xmin = 3.95, 
           xmax = 4.02, 
           ymin = 0,
           ymax = Inf,
           fill = custom_red, colour = custom_red, alpha = 0.1) +
  annotate(geom = "text",
           x = 3.75,
           y = 210,
           label = "BB B") +
  annotate(geom = "rect", 
           xmin = 3.2, 
           xmax = 3.3, 
           ymin = 0,
           ymax = Inf,
           fill = custom_red, colour = custom_red, alpha = 0.1) +
  annotate(geom = "text",
           x = 3,
           y = 210,
           label = "BB C") +
  geom_area(alpha = 0.3, 
            position = "identity") +
  facet_grid(sample ~ .,
             labeller = labeller(sample = c("kb_liquid" = "KB liquid",
                                            "kb_solid" = "KB solid",
                                            "lb_liquid" = "LB liquid",
                                            "lb_solid" = "LB solid"))) +
  scale_x_continuous(limits = c(0,10), breaks = c(0:10)) +
  scale_fill_manual(values = c(custom_red, brew_blue, custom_grey),
                    labels = c("Replicate 1", "Replicate 2", "Replicate 3"),
                    name = "") +
  scale_color_manual(values = c(custom_red, brew_blue, custom_grey),
                     labels = c("Replicate 1", "Replicate 2", "Replicate 3"),
                     name = "") +
  theme_bw() +
  labs(y = "A312 (mAU)",
       x = "Retention time (min)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "bottom",
    plot.margin = margin(t = 5)
  )

# AUC faceted on compound
auc_p_compound <-  auc %>%
  dplyr::filter(str_detect(compound, "bb"),
                str_detect(sample, "solid")) %>%
  mutate(compound = as_factor(compound),
         compound = fct_relevel(compound, c("bb_c", "bb_b", "bb"))) %>% 
  ggplot(aes(sample, auc_uv_312+1)) +
  stat_summary(fun.data = mean_se ,
               geom = "errorbar",
               col = custom_red,
               width = 0.2,
               linewidth = 1) +
  geom_point(alpha = 0.5) +
  geom_pwc(method = "t_test",
           y.position = 3.3) +
  facet_grid(. ~ compound,
             scales = "free",
             labeller = labeller(compound = c("bb" = "Bacillibactin\n RT = 4.25 - 4.45 mins",
                                              "bb_b" = "Bacillibactin B\n RT = 3.99 mins",
                                              "bb_c" = "Bacillibactin C\n RT = 3.25 mins"))) +
  scale_y_log10(limits = c(1, 5*10^3)) +
  scale_x_discrete(labels = c("KB solid", "LB solid", "KB liquid", "LB liquid")) +
  theme_bw() +
  labs(y = "AUC (A312)",
       x = "") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "bottom",
    plot.margin = margin(t = 0,  # Top margin
                         r = 0,  # Right margin
                         b = 0,  # Bottom margin
                         l = 0) # Left margin
  )
  
# LCMS with inset
dad_p_inset <- dad_p +
  inset_element(auc_p_compound,
                left = 0.6,
                bottom = 0.04,
                right = 0.99,
                top = 0.49)

# Figure 1 plots
(ratio_p + cfu_p) / dad_p_inset

ggsave(filename = "fig1_plots.pdf",
       path = here("figures", "figure_1_antagonism"),
       dpi = 600,
       scale = 2,
       width = 185,
       height = 140,
       units = "mm"
)

# Presentation
cfu_p <- cfu_p +
  theme(text = element_text(size = 15),
        strip.text = element_text(size = 15))

ratio_p <- ratio_p +
  theme(text = element_text(size = 15),
        strip.text = element_text(size = 15))


plots_aligned <- align_patches(cfu_p, ratio_p)
plots_aligned[[2]] 
ggsave(filename = "fig1b.png",
       path = here("figures", "presentation"),
       dpi = 600,
       scale = 2,
       width = 185,
       height = 100,
       units = "mm")

plots_aligned[[1]] 
ggsave(filename = "fig1c.png",
       path = here("figures", "presentation"),
       dpi = 600,
       scale = 2,
       width = 185,
       height = 100,
       units = "mm")



dad_p <- dad %>% 
  ggplot(aes(x = time_m, 
             y = value, 
             fill = fct_reorder(rep, value, .desc = F),
             col = fct_reorder(rep, value, .desc = F))) +
  annotate(geom = "rect", 
           xmin = 4.1, 
           xmax = 4.5, 
           ymin = 0,
           ymax = Inf,
           fill = custom_red, colour = custom_red, alpha = 0.1) + 
  annotate(geom = "text",
           x = 4.65,
           y = 210,
           label = "BB") +
  annotate(geom = "rect", 
           xmin = 3.95, 
           xmax = 4.02, 
           ymin = 0,
           ymax = Inf,
           fill = custom_red, colour = custom_red, alpha = 0.1) +
  annotate(geom = "text",
           x = 3.75,
           y = 210,
           label = "BB B") +
  annotate(geom = "rect", 
           xmin = 3.2, 
           xmax = 3.3, 
           ymin = 0,
           ymax = Inf,
           fill = custom_red, colour = custom_red, alpha = 0.1) +
  annotate(geom = "text",
           x = 3,
           y = 210,
           label = "BB C") +
  geom_area(alpha = 0.3, 
            position = "identity") +
  facet_grid(sample ~ .,
             labeller = labeller(sample = c("kb_liquid" = "KB liquid",
                                            "kb_solid" = "KB solid",
                                            "lb_liquid" = "LB liquid",
                                            "lb_solid" = "LB solid"))) +
  scale_x_continuous(limits = c(0,10), breaks = c(0:10)) +
  scale_fill_manual(values = c(custom_red, brew_blue, custom_grey),
                    labels = c("Replicate 1", "Replicate 2", "Replicate 3"),
                    name = "") +
  scale_color_manual(values = c(custom_red, brew_blue, custom_grey),
                     labels = c("Replicate 1", "Replicate 2", "Replicate 3"),
                     name = "") +
  theme_bw() +
  labs(y = "A312 (mAU)",
       x = "Retention time (min)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    text = element_text(size = 15),
    legend.position = "bottom",
    plot.margin = margin(t = 5)
  )

dad_p
ggsave(filename = "fig1d.png",
       path = here("figures", "presentation"),
       dpi = 600,
       scale = 1.5,
       width = 210,
       height = 150,
       units = "mm")

library(multcompView)
auc_aov <- auc %>%
  dplyr::filter(str_detect(compound, "bb")) %>%
  group_by(compound) %>% 
  nest() %>% 
  mutate(model = map(data, ~aov(auc_uv_312 ~ sample, data = .))) %>%
  mutate(tukey = map(model, TukeyHSD)) %>%  
  dplyr::select(-data)

auc_letters <- auc_aov %>%
  mutate(letters = map(tukey, function(x) {
    multcompLetters(x$sample[, "p adj"], reversed = T)$Letters
    })
    ) %>% 
  mutate(letters = map(letters, enframe)) %>% 
  unnest(letters) %>% 
  rename(sample = name)


# AUC faceted on compound
auc_p_compound <-  auc %>%
  dplyr::filter(str_detect(compound, "bb")) %>%
  mutate(compound = as_factor(compound),
         compound = fct_relevel(compound, c("bb_c", "bb_b", "bb"))) %>% 
  ggplot(aes(sample, auc_uv_312+1)) +
  stat_summary(fun.data = mean_se ,
               geom = "errorbar",
               col = custom_red,
               width = 0.2,
               linewidth = 1) +
  geom_point(alpha = 0.5) +
  geom_text(data = auc_letters,
            aes(x = sample,
                label = value),
            y = 3.5,
            size = 5,
            inherit.aes = F) +
  facet_grid(. ~ compound,
             scales = "free",
             labeller = labeller(compound = c("bb" = "Bacillibactin\n RT = 4.25 - 4.45 mins",
                                              "bb_b" = "Bacillibactin B\n RT = 3.99 mins",
                                              "bb_c" = "Bacillibactin C\n RT = 3.25 mins"))) +
  scale_y_log10(limits = c(1, 5*10^3)) +
  scale_x_discrete(labels = c("KB solid", "LB solid", "KB liquid", "LB liquid")) +
  theme_bw() +
  labs(y = "AUC (A312)",
       x = "") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    text = element_text(size = 15),
    legend.position = "bottom",
    plot.margin = margin(t = 0,  # Top margin
                         r = 0,  # Right margin
                         b = 0,  # Bottom margin
                         l = 0) # Left margin
  )

auc_p_compound
ggsave(filename = "fig1d_inset.png",
       path = here("figures", "presentation"),
       dpi = 600,
       scale = 2.5,
       width = 120,
       height = 60,
       units = "mm")




#### Figure S1: Antagonism is medium-specific and not due to space restraints ####
# Data import
  # Growth curves
data <- map(list.files(here("data", "growth curves"), pattern = ".xlsx", full.names = T), function(x) readxl::read_xlsx(x, col_names = F))
md_list <- map(list.files(here("data", "growth curves"), pattern = ".csv", full.names = T), read_csv2)

data_parsed <- map2(data, md_list, function(x, md) {
  
  data_idx <- which(grepl("Time", pull(x[,2])))
  
  df_list <- map(data_idx, function(idx) {
    block_end <- idx + min(
      which(
        pull(
          x[idx:nrow(x), 4]) 
        %in% c("", NA)
      )
    )
    
    df <- x %>%
      slice(idx:block_end) %>%
      dplyr::select(-1:-3)
    
    names(df) <- as.character(df[1, ])
    df <- df %>% 
      slice(-1) %>% 
      pivot_longer(cols = everything(), 
                   names_to = "well") %>% 
      arrange(well) %>% 
      group_by(well) %>% 
      mutate(value = as.numeric(value),
             measure = pull(x[idx-2, 1]) %>% 
               str_remove(":.*"),
             time_m = seq(from = 0,
                          to = (pluck(tally(.), 2, 1)-1) * 10,
                          by = 10),
             time_h = time_m/60) %>% 
      ungroup()
  })
  
  df <- do.call("rbind", df_list) %>% 
    left_join(md, by = "well") %>% 
    dplyr::filter(!is.na(value),
           sample != "cont")
  
})

df_growth <- do.call("rbind", data_parsed)

# Subtract blanks (the minimum values of each well)
df_blk <- df_growth %>%
  group_by(well, measure) %>% 
  mutate(value_blk = value - min(value),
         measure = factor(measure, levels = c("od", "rfp", "gfp"))) %>% 
  dplyr::filter(!sample == "blk")

# Plot
# Growth curves
gc_p <- df_blk %>% 
  ggplot(aes(x = time_h, y = value_blk, col = sample, linetype = sample)) +
  geom_smooth(linewidth = 0.7) +
  facet_grid(measure ~ medium,
             scales = "free",
             labeller = labeller(measure = c(od = "OD600", rfp = "RFP (RFU)", gfp = "GFP (RFU)"),
                                 medium = c(kb = "King's B", lb = "LB"))) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_brewer(palette = "Paired",
                     name = "Culture",
                     labels = c("3610",
                                "3610 + PS92",
                                "dhbA",
                                "dhbA + PS92",
                                "PS92")) +
  scale_linetype_manual(name = "Culture",
                        values = c(2,1,2,1,2),
                        labels = c("3610",
                                   "3610 + PS92",
                                   "dhbA",
                                   "dhbA + PS92",
                                   "PS92")
  ) +
  theme_bw() +
  labs(y = NULL,
       x = "Time (h)") +
  theme(text = element_text(size = text_size),
        strip.background = element_blank(),
        strip.text = element_text(size = text_size)) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

# Colony area ratio
area_p <- df %>% 
  ggplot(aes(x = strain, y = area_mm2, fill = type)) +
  geom_violin() +
  ggforce::geom_sina(alpha = 0.5, size = 0.5) +
  stat_compare_means(inherit.aes = T,
                     method = "t.test",
                     label = "p.format",
                     paired = F) +
  facet_grid(interaction ~ time,
             labeller = labeller(interaction = label_parsed,
                                 time = label_value)) +
  scale_y_continuous() +
  scale_fill_manual(name = "",
                    labels = c("Coculture", "Monoculture"),
                    values = c(brew_blue, custom_red)) +
  scale_x_discrete(labels = c("Bacillus",
                              "Pseudomonas"
                              )
  ) +
  theme_bw() +
  labs(y = expression(paste("Area (", mm^2, ")")),
       x = "") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size)
  )

# Combine
area_p + gc_p

ggsave(filename = "figS1_plots.png",
       path = here("figures", "base_interaction"),
       dpi = 600,
       scale = 2,
       height = 100,
       width = 185,
       units = "mm"
)

# Presentation
gc_p <- df_blk %>% 
  dplyr::filter(medium == "kb") %>% 
  ggplot(aes(x = time_h, y = value_blk, col = sample, linetype = sample)) +
  geom_smooth(linewidth = 0.7) +
  facet_wrap(medium ~ measure,
             scales = "free_y",
             labeller = labeller(measure = c(od = "OD600", rfp = "RFP (RFU)", gfp = "GFP (RFU)"),
                                 medium = c(kb = "King's B", lb = "LB"))) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_brewer(palette = "Paired",
                     name = "Culture",
                     labels = c("3610",
                                "3610 + PS92",
                                "dhbA",
                                "dhbA + PS92",
                                "PS92")) +
  scale_linetype_manual(name = "Culture",
                        values = c(2,1,2,1,2),
                        labels = c("3610",
                                   "3610 + PS92",
                                   "dhbA",
                                   "dhbA + PS92",
                                   "PS92")
  ) +
  theme_bw() +
  labs(y = NULL,
       x = "Time (h)") +
  theme(text = element_text(size = 15),
        strip.background = element_blank(),
        strip.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

gc_p
ggsave(filename = "figS1b.png",
       path = here("figures", "presentation"),
       dpi = 600,
       scale = 2,
       width = 185,
       height = 100,
       units = "mm")

#### Figure 2: BB intermediates are public goods ####
# Data import
  # dhb mutants
df_dhb <- read.csv(here("data", "area", "dhb_mutants", "dhb_mutants_area.csv"), sep = ",") %>% 
  mutate(Label = str_remove(Label, "-[:digit:]{2}"),
         Label = str_remove(Label, ".czi")) %>% 
  separate(Label, into = c("interaction", "type", "strain")) %>% 
  dplyr::select(-X) %>% 
  rename(area_um2 = Area) %>% 
  mutate(strain = ifelse(type == "co", 
                         strain,
                         ifelse(type == "bs", "bs", "ps")),
         type = ifelse(type == "co", type, "mono"),
         area_mm2 = area_um2/10^6,
         var = str_c(strain, type, sep = "_"),
         interaction = as.factor(interaction))

split_list <- df_dhb %>% 
  group_by(type) %>% 
  group_split()

ratio_dhb <- split_list[[2]] %>%
  group_by(strain, interaction) %>% 
  summarize(mono_mean = mean(area_mm2)) %>% 
  right_join(ungroup(split_list[[1]]), by = c("strain", "interaction")) %>% 
  mutate(l2fc = log2(area_mm2 / mono_mean),
         ratio = area_mm2/mono_mean,
         strain = as.factor(strain))


dhb_test <- compare_means(area_mm2 ~ type, 
                          df_dhb, 
                          method = "t.test", 
                          group.by = c("strain", "interaction")) %>% 
  mutate(p.signif = ifelse(p.adj > 0.05, "ns", p.signif),
         strain = as.factor(strain))


levels(ratio_dhb$strain) = c("bs" = expression(paste(italic(B.~subtilis))),
                             "ps" = expression(paste(italic(P.~marginalis)))
)
levels(dhb_test$strain) = c("bs" = expression(paste(italic(B.~subtilis))),
                            "ps" = expression(paste(italic(P.~marginalis)))
)

  # complemented dhb mutants
df_com <- read.csv(here("data", "area", "complementation_area.csv"), sep = ",") %>% 
  mutate(Label = str_remove(Label, "-[:digit:]{2}"),
         Label = str_remove(Label, ".czi")) %>% 
  separate(Label, into = c("interaction", "inoc_ratio", "type", "strain")) %>% 
  dplyr::select(-X) %>% 
  rename(area_um2 = Area) %>% 
  mutate(strain = ifelse(type == "co", 
                         strain,
                         ifelse(type == "bs", "bs", "ps")),
         type = ifelse(type == "co", type, "mono"),
         area_mm2 = area_um2/10^6,
         var = str_c(strain, type, sep = "_"),
         interaction = as.factor(interaction))

df_wtcom <- df_com %>% 
  dplyr::filter(str_detect(interaction, "bs"))

df_mtcom <- df_com %>% 
  dplyr::filter(str_detect(interaction, "bs", negate = T)) %>% 
  mutate(strain = as_factor(strain))

levels(df_wtcom$interaction) = c("bsdhbA" = expression(paste(" WT", 
                                                             " and ", 
                                                             ~Delta, italic(dhbA))),
                                 "bsdhbF" = expression(paste(" WT", 
                                                             " and ", 
                                                             ~Delta, italic(dhbF))),
                                 "dhbadhbb" = "dhbadhbb",
                                 "dhbadhbf" = "dhbadhbf",
                                 "dhbedhbf" = "dhbedhbf"
)

levels(df_mtcom$strain) = c("bs" = expression(paste(italic(B.~subtilis))),
                            "ps" = expression(paste(italic(P.~marginalis)))
)

split_list <- df_mtcom %>% 
  group_by(type) %>% 
  group_split()

ratio_mtcom <- split_list[[2]] %>%
  group_by(interaction, strain) %>% 
  summarize(mono_mean = mean(area_mm2), .groups = "drop") %>% 
  right_join(ungroup(split_list[[1]]), by = c("interaction", "strain")) %>% 
  mutate(l2fc = log2(area_mm2 / mono_mean),
         ratio = area_mm2/mono_mean)


compl_test <- compare_means(area_mm2 ~ type, 
                            df_mtcom, 
                            method = "t.test", 
                            group.by = c("strain", "interaction")) %>% 
  mutate(p.signif = ifelse(p.adj > 0.05, "ns", p.signif),
         strain = as.factor(strain))

  
# Plot
# dhb mutant interactions
dhb_p <- ratio_dhb %>%  
  ratio_plot_fill(aes(interaction, ratio, fill = strain), 
                  text_pos = 4) +
  # geom_text(data = dhb_test, 
  #           aes(x = interaction, label = p.signif),
  #           y = 1.2,
  #           size = 5) +
  scale_x_discrete(labels = c("dhbB",
                              "dhbC",
                              "dhbE",
                              "dhbF")) +
  scale_fill_manual(values = c(cust_mag, cust_green),
                    guide = "none") +
  facet_grid(.~strain,
             labeller = labeller(strain = label_parsed)) +
  labs(y = "Area coculture / area monoculture",
       x = "")

# Complemented dhb mutants
compl_p <- ratio_mtcom %>% 
  ratio_plot_fill(aes(x = interaction, 
                      y = ratio,
                      fill = strain), 
                  text_pos = 3) +
  # geom_text(data = compl_test,
  #           aes(x = interaction, label = p.signif),
  #           y = 1.2,
  #           size = 5) +
  facet_grid(. ~ strain,
             labeller = label_parsed) +
  scale_fill_manual(values = c(cust_mag, cust_green),
                    guide = "none") +
  scale_x_discrete(name = "",
                   labels = c(expression(paste(italic(dhbA), " + ", italic(dhbB))),
                              expression(paste(italic(dhbA), " + ", italic(dhbF))),
                              expression(paste(italic(dhbE), " + ", italic(dhbF)))
                   )
  )

# Combined
(dhb_p + compl_p) +
  plot_annotation(tag_levels = 'A')

ggsave(filename = "fig2_plots.png",
       path = here("figures", "figure_2_publicgoods"),
       dpi = 600,
       scale = 2,
       width = 185,
       height = 60,
       units = "mm"
)

# Presentation
dhb_p + compl_p &
  theme(text = element_text(size = 15),
        strip.text = element_text(size = 15))

ggsave(filename = "fig2ab.png",
       path = here("figures", "presentation"),
       dpi = 600,
       scale = 2.5,
       width = 200,
       height = 100,
       units = "mm"
)

#### Figure S2: dhb mutants can be complemented by 3610 WT ####

#### Figure 3: Pseudomonas size restriction is related to sfp and feuA ####
# Data import
  # Mutant interaction colony areas
df_mutant <- map(list.files(here("data", "area", "mutants"), "area.csv", full.names = T), function(x) read.csv(x, sep = ";")) %>% 
  map2(list.files(here("data", "area", "mutants"), "area.csv"), function(x, y) {
    x %>%
      mutate(run = str_remove(y, "_area.csv"),
             Label = str_remove(Label, "-[:digit:]{2}"),
             Label = str_remove(Label, ".czi")) %>% 
      separate(Label, into = c("interaction", "type", "strain")) %>% 
      dplyr::select(-X) %>% 
      rename(area_um2 = Area) %>% 
      mutate(strain = ifelse(type == "co", 
                             strain,
                             ifelse(type == "bs", "bs", "ps")),
             type = ifelse(type == "co", type, "mono"),
             area_mm2 = area_um2/10^6,
             var = str_c(strain, type, sep = "_"),
             interaction = as.factor(interaction),
      ) %>% 
      dplyr::select(-strain, -var)
  })

# All monocultures are Pseudomonas (I didn't measure Bs)
df_mutant <- do.call("rbind", df_mutant) %>% 
  mutate(interaction = ifelse(type == "mono", "ps", as.character(interaction)),
         interaction = ifelse(interaction == "ps", interaction, str_remove(interaction, "ps")),
  ) 

# Split by mono/co-culture and calculate ratio
split_list <- df_mutant %>% 
  group_by(type) %>% 
  group_split()

ratio_mutant <- split_list[[2]] %>% 
  group_by(run) %>% 
  summarize(mono_mean = mean(area_mm2)) %>% 
  right_join(ungroup(split_list[[1]]), by = c("run")) %>% 
  mutate(l2fc = log2(area_mm2 / mono_mean),
         ratio = area_mm2/mono_mean,
         interaction = str_c(str_extract(interaction, "[:alpha:]{3}"), 
                             str_to_upper(str_extract(interaction, pattern = "(?<=[:alpha:]{3}).*")))) %>% 
  left_join(read_csv2(here("data/area/mutants", "md.csv")), by = "interaction")


# non-parametric test comparing coculture size to monoculture size from only the mutant experiments
# Does not test on ratio-transformed data
mutant_test <- compare_means(area_mm2 ~ interaction, df_mutant, method = "wilcox.test", ref.group = "ps", group.by = c("run")) %>% 
  mutate(p.signif = ifelse(p.adj > 0.05, "ns", p.signif),
         group2 = as_factor(group2)) %>% 
  rename("interaction" = group2) %>% 
  # Complicated way to make last letter in gene uppercase
  mutate(interaction = str_c(str_extract(interaction, "[:alpha:]{3}"), 
                             str_to_upper(str_extract(interaction, pattern = "(?<=[:alpha:]{3}).*"))),
         interaction = as_factor(interaction),
         interaction = fct_relevel(interaction, 
                                   levels(fct_reorder(ratio_mutant$interaction, ratio_mutant$ratio, .desc = T))))

# PS92 colony size with mutants
ratio_mutant %>% 
  ggplot(aes(fct_reorder(interaction, ratio, .desc = T), 
             y = ratio,
             fill = BB_effect)) +
  geom_hline(data = dplyr::filter(reference, str_detect(var, "ps"), interaction == "bsps"),
             aes(yintercept = mean),
             linetype = 2) +
  geom_rect(data = dplyr::filter(reference, str_detect(var, "ps"), interaction == "bsps"),
            aes(ymax = max,
                ymin = min),
            xmin = -Inf, 
            xmax = Inf,
            alpha = .1,
            fill = "black",
            inherit.aes = F) +
  geom_hline(data = dplyr::filter(reference, str_detect(var, "ps"), interaction == "dhbps"),
             aes(yintercept = mean),
             linetype = 2,
             col = custom_red) +
  geom_rect(data = dplyr::filter(reference, str_detect(var, "ps"), interaction == "dhbps"),
            aes(ymax = max,
                ymin = min),
            xmin = -Inf, 
            xmax = Inf,
            alpha = .1,
            fill = custom_red,
            inherit.aes = F) +
  geom_text(data = dplyr::filter(reference, str_detect(var, "ps"), interaction == "bsps"),
            aes(y = ifelse(max < 1, 
                           max - 0.3, 
                           max - 0.04)),
            x = 21 + 0.5,
            label = "Bs vs Ps",
            hjust = 1,
            size = 5,
            fontface = "bold",
            inherit.aes = F) +
  geom_text(data = dplyr::filter(reference, str_detect(var, "ps"), interaction == "dhbps"),
            aes(y = max - 0.04),
            x = 21 + 0.5,
            label = "dhbA vs Ps",
            hjust = 1,
            size = 5,
            fontface = "bold",
            col = custom_red,
            inherit.aes = F) +
  geom_boxplot(outlier.shape = NA) +
  ggforce::geom_sina() +
  scale_y_continuous(limits = c(0.1, 1.7),
                     breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6)) +
  scale_fill_manual(values = c(custom_red, brew_blue, custom_grey),
                    labels = c("Less BB", "More BB/no difference", "Other/Unknown"),
                    name = "Mutant outcome") +
  theme_bw() +
  labs(y = "Area coculture / area monoculture") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "bottom"
  ) +
  geom_text(data = mutant_test, 
            aes(x = interaction, 
                label = p.signif),
            y = 1.2,
            size = text_size,
            inherit.aes = F) +  
  labs(x = "Mutant partner")

ggsave(filename = "mutant_ratio_labels.png",
       device = png,
       path = here("figures"),
       dpi = 600,
       width = 13,
       height = 8
)


# Presentation

#### Figure 4: Bacillus offense and defense are iron-dependent ####
# Data import
df_fe <- map(list.files(here("data", "area", "iron_exp"), ".csv", full.names = T), function(x) read.csv(x, sep = ";")) %>% 
  map2(list.files(here("data", "area", "iron_exp"), ".csv"), function(x, y) {
    x %>%
      mutate(run = str_remove(y, "_area.csv"),
             Label = str_remove(Label, "-[:digit:]{2}"),
             Label = str_remove(Label, ".czi")) %>% 
      separate(Label, into = c("media", "conc", "interaction", "type", "strain"), sep = "_") %>% 
      dplyr::select(-X) %>% 
      rename(area_um2 = Area) %>% 
      mutate(strain = ifelse(type == "co", 
                             strain,
                             ifelse(type == "bs", "bs", "ps")),
             type = ifelse(type == "co", type, "mono"),
             area_mm2 = area_um2/10^6,
             var = str_c(strain, type, sep = "_"),
             interaction = as.factor(interaction))
  })

df_fe <- do.call("rbind", df_fe)

# Split data by mono/co-culture and calculate ratio
split_list <- df_fe %>% 
  group_by(type) %>% 
  group_split()

ratio_fe <- split_list[[2]] %>%
  group_by(strain, conc, media, run, interaction) %>% 
  summarize(mono_mean = mean(area_mm2), .groups = "drop") %>% 
  right_join(ungroup(split_list[[1]]), by = c("strain", "interaction", "conc", "media", "run")) %>% 
  mutate(l2fc = log2(area_mm2 / mono_mean),
         ratio = area_mm2/mono_mean,
         conc = ifelse(media == "fecit",
                       str_c(conc, " µg/mL"),
                       str_c(conc, " µM")),
         conc = forcats::fct_relevel(conc, c("20 µg/mL", "25 µM", "50 µM", "100 µg/mL", "100 µM", "500 µg/mL")),
         strain = as.factor(strain))

ratio_fe_labels = c("kb-bpd" = "King's B with BPD",
                    "lb-bpd" = "LB with BPD",
                    "fecit" = "King's B with ferric citrate",
                    "fecl" = "King's B with ferric chloride"
)

levels(ratio_fe$strain) = c("bs" = expression(paste(italic(B.~subtilis))),
                            "ps" = expression(paste(italic(P.~marginalis)))
)

# Plot
# BPD addition
bpd_p <- ratio_fe %>% 
  dplyr::filter(str_detect(media, "kb-bpd")) %>% 
  ratio_plot_gradient(aes(interaction, ratio, fill = conc), 
                      text_pos = 2) +
  facet_grid(. ~ strain,
             labeller = label_parsed) +
  scale_fill_brewer(palette = "Blues",
                    name = "BPD\nConcentration") +
  labs(x = "",
       title = "King's B +BPD") +
  scale_x_discrete(labels = c("bsps" = "Wildtype",
                              "dhbps" = "dhbA",
                              "besaps" = "besA",
                              "feuaps" = "feuA"))

# FeCit addition
fecit_p <- ratio_fe %>% 
  dplyr::filter(str_detect(media, "fecit")) %>% 
  ratio_plot_gradient(aes(interaction, ratio, fill = conc), 
                      text_pos = 4) +
  facet_grid(. ~ strain,
             labeller = label_parsed) +
  scale_fill_brewer(palette = "Blues",
                    name = "FeCit\nConcentration") +
  labs(x = "",
       title = "King's B +Ferric citrate",
       y = "") +
  scale_x_discrete(labels = c("bsps" = "Wildtype",
                              "dhbps" = "dhbA",
                              "besaps" = "besA",
                              "feuaps" = "feuA"))

# FeCl addition
fecl_p <- ratio_fe %>% 
  dplyr::filter(str_detect(media, "fecl")) %>% 
  ratio_plot_gradient(aes(interaction, ratio, fill = conc), 
                      text_pos = 4) +
  facet_grid(. ~ strain,
             labeller = label_parsed) +
  scale_fill_brewer(palette = "Blues",
                    name = "FeCl\nConcentration") +
  labs(x = "Bacillus partner",
       title = "King's B +Ferric chloride",
       y = "") +
  scale_x_discrete(labels = c("bsps" = "Wildtype",
                              "dhbps" = "dhbA",
                              "besaps" = "besA",
                              "feuaps" = "feuA"))

# Collect vertically
bpd_p / fecit_p / fecl_p + 
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = margin(r = 0, l = 0))

ggsave(filename = "fig3_plots.pdf",
       path = here("figures", "figure_3_irondepend"), 
       scale = 2, 
       dpi = 600,
       width = 80,
       height = 180,
       unit = "mm")

# Presentation
ratio_fe %>% 
  dplyr::filter(str_detect(media, "kb-bpd")) %>% 
  ratio_plot_gradient(aes(interaction, ratio, fill = conc), 
                      text_pos = 2) +
  facet_grid(. ~ strain,
             labeller = label_parsed) +
  scale_fill_brewer(palette = "Blues",
                    name = "BPD\nConcentration") +
  labs(x = "",
       title = "King's B -Fe(II)") +
  scale_x_discrete(labels = c("bsps" = "Wildtype",
                              "dhbps" = "dhbA",
                              "besaps" = "besA",
                              "feuaps" = "feuA")) +
  theme(text = element_text(size = 15),
        strip.text = element_text(size = 15))

ggsave(filename = "fig3a.png",
       path = here("figures", "presentation"), 
       scale = 2, 
       dpi = 600,
       width = 120,
       height = 90,
       unit = "mm")

fecit_p + fecl_p &
  theme(text = element_text(size = 15),
        strip.text = element_text(size = 15))  

ggsave(filename = "fig3bc.png",
       path = here("figures", "presentation"), 
       scale = 2, 
       dpi = 600,
       width = 180,
       height = 90,
       unit = "mm")



#### Figure S3: Media components influence iron dependency and antagonism ####

#### Figure 5: dhbA is unable to affect PS92 transcriptome ####
# Data import
  # Kallisto counts
files <- list.files(here("data", "rnaseq_raw", "bakta"), "*.h5", recursive = T, full.names = T) 
names(files) <- list.files(here("data", "rnaseq_raw", "bakta"), recursive = F)

txi_kallisto <- tximport(files, type = "kallisto", txOut = TRUE, tx2gene = T) %>%
  map(function(x) {
    if(is.matrix(x)) {
      return(x[!str_detect(row.names(x), "005363|001014"), ])
    } else return(x)
  })

txi_bs <- map(txi_kallisto, function(x) {
  if(is.matrix(x)) {
    return(x[str_detect(row.names(x), "B4U62"), -13:-15])
  } else return(x)
})

txi_ps <- map(txi_kallisto, function(x) {
  if(is.matrix(x)) {
    return(x[str_detect(row.names(x), "FNPKGJ"), c(-1:-3, -7:-9)])
  } else return(x)
})

# Features from BAKTA-annotated GFFs
feats <- rbind(rtracklayer::import.gff3(here("md", "PS92_bakta.gff3")) %>%
                 as.data.frame() %>% 
                 dplyr::filter(type %in% c("CDS", "ncRNA", "regulatory_region")) %>% 
                 dplyr::select(seqnames, type, ID, Name, gene, locus_tag, product) %>% 
                 mutate(Name = str_remove(Name, "extdb:")),
               rtracklayer::import.gff3(here("md", "3610_comb.gff3")) %>%
                 as_tibble() %>% 
                 dplyr::filter(type %in% c("CDS", "ncRNA", "regulatory_region")) %>% 
                 dplyr::select(seqnames, type, ID, Name, gene, locus_tag, product)
)

# Sample metadata
md <- read_csv2(here("md", "rnaseq_md.csv"))
rownames(md) <- colnames(txi_kallisto$counts)

md_bs <- md %>%
  dplyr::filter(interaction != "ps")

md_ps <- md %>%
  dplyr::filter(str_detect(interaction, "ps"))

#### DESeq2
dds_bs <- DESeqDataSetFromTximport(txi_bs, colData = md_bs, design = ~interaction)
dds_bs <- dds_bs[rowSums(counts(dds_bs)) >= 10, ] %>%
  DESeq()

dds_ps <- DESeqDataSetFromTximport(txi_ps, colData = md_ps, design = ~interaction)
dds_ps <- dds_ps[rowSums(counts(dds_ps)) >= 10, ] %>%
  DESeq()

# High l2fc means up in numerator, low l2fc means up in denominator
# Co in numerator, mono in denominator: high l2fc means up in coculture
comps <- list(bs_covsmono = results(dds_bs, contrast = c("interaction", "bs_ps", "bs"), 
                                    tidy = T, pAdjustMethod = p_adj_meth),
              bs_mtvswt = results(dds_bs, contrast = c("interaction", "mt", "bs"), 
                                  tidy = T, pAdjustMethod = p_adj_meth),
              mt_covsmono = results(dds_bs, contrast = c("interaction", "mt_ps", "mt"), 
                                    tidy = T, pAdjustMethod = p_adj_meth),
              bs_co_mtvswt = results(dds_bs, contrast = c("interaction", "mt_ps", "bs_ps"), 
                                     tidy = T, pAdjustMethod = p_adj_meth),
              ps_covsmono_wt = results(dds_ps, contrast = c("interaction", "bs_ps", "ps"),
                                       tidy = T, pAdjustMethod = p_adj_meth),
              ps_covsmono_mt =  results(dds_ps, contrast = c("interaction", "mt_ps", "ps"),
                                        tidy = T, pAdjustMethod = p_adj_meth),
              ps_co_mtvswt = results(dds_ps, contrast = c("interaction", "mt_ps", "bs_ps"),
                                     tidy = T, pAdjustMethod = p_adj_meth)
)

# Clean-up
rm(list = ls(pattern = "txi"))

# Data wrangle and custom annotation
comps <- map(comps, function(df) {
  df %>% 
    as_tibble() %>% 
    mutate(p_label = ifelse(padj < p_thresh, "yes", "no"),
           lfc_label = ifelse(abs(log2FoldChange) > lfc_thresh, "yes", "no"),
           row = ifelse(grepl("B4U62", row), str_extract(row, pattern = "B4U62_[:digit:]*"), row),
           row = ifelse(grepl("FNPKGJ", row), str_extract(row, pattern = "FNPKGJ_[:digit:]*"), row)) %>% 
    dplyr::rename(locus_tag = row) %>% 
    left_join(., feats, by = "locus_tag") %>% 
    mutate(change = ifelse(log2FoldChange > 0, "up", "down"))
})


#### Functional annotation
# Append functional annotation from Eggnog-mapper
ps_egg <- read_tsv(here("md", "ps92_emapper.tsv"), skip = 4) %>% 
  dplyr::rename(locus_tag = query, gene_emap = Preferred_name) %>% 
  mutate(locus_tag = str_remove(locus_tag, "[:alpha:].*-"))

bs_egg <- read_tsv(here("md", "3610_emapper.tsv"), skip = 4) %>% 
  dplyr::rename(locus_tag = query, gene_emap = Preferred_name) %>% 
  mutate(locus_tag = str_remove(locus_tag, "[:alpha:].*-"))

comps <- map(comps, function(comp) {
  if (any(grepl("B4U62", comp$locus_tag))) {
    comp %>% 
      left_join(bs_egg, by = "locus_tag") %>% 
      mutate(gene = if_else(gene == "ybgJ", "glsA", gene), # Manual http://subtiwiki.uni-goettingen.de/v4/gene?id=1A415DC85354373EA4866733C3AE43F510BD3C56
             gene = if_else(gene == "gltD", "gltB", gene), # Manual (gltD does not exist in Subtiwiki)
             gene = if_else(gene == "yolJ", "sunS", gene)) # Manual http://subtiwiki.uni-goettingen.de/v4/gene?id=6B2DC9873B7B693D9D0FC0D773740140D0640C9E
  } else if (any(grepl("FNPKGJ", comp$locus_tag))) {
    comp %>% 
      left_join(ps_egg, by = "locus_tag")
  }
})

# Pivot data - let each contrast be a variable and each gene an observation
comps_wide <- map2(comps, names(comps), function(df, contrast) {
  l2fc_name = str_c(contrast, "_l2fc")
  padj_name = str_c(contrast, "_padj")
  change_name = str_c(contrast, "_change")
  df %>%
    dplyr::select(locus_tag, log2FoldChange, padj, Name, product, gene, gene_emap, change) %>% 
    dplyr::rename(!!l2fc_name := log2FoldChange,
                  !!padj_name := padj,
                  !!change_name := change) %>% 
    mutate(gene_emap = ifelse(gene_emap == "-", NA, gene_emap))
})

comps_wide_ps <- comps_wide %>% 
  purrr::keep_at(str_detect(names(comps), "ps")) %>% 
  map(function(df) {
    df %>% # Manual annotation
      mutate(gene_label = ifelse(str_detect(Name, "[Aa]lg"), "Alginate biosynthesis", NA),
             gene_label = ifelse(str_detect(gene_emap, "pvd"), "Pyoverdine biosynthesis", gene_label),
             gene_label = ifelse(str_detect(Name, "type VI"), "Type VI secretion system", gene_label),
             gene_label = ifelse(str_detect(Name, ".*[Gg]ac.*|.*[Rr]sm.*"), "Gac/Rsm", gene_label))
  }) %>% 
  purrr::reduce(full_join, 
                by = c("locus_tag", "Name", "gene", "gene_emap", "gene_label"))

comps_wide_bs <- comps_wide[1:4] %>% 
  purrr::reduce(full_join, 
                by = c("locus_tag", "Name", "product", "gene", "gene_emap")) %>% 
  # Manual annotation
  mutate(gene_label = ifelse(str_detect(product, regex("^spo| spo")), "Sporulation", NA),
         gene_label = ifelse(str_detect(product, "ger"), "Germination", gene_label),
         gene_label = ifelse(str_detect(gene, "dhb"), "Bacillibactin", gene_label),
         gene_label = ifelse(str_detect(gene, "eps|tas|tap|sip|rap|bsl"), "Biofilm", gene_label),
         gene_label = ifelse(str_detect(gene, "srf|sbo|pks|pps.|sun|^bac|sqh|^cyp|^yvm"), "Secondary metabolites", gene_label),
         gene_label = ifelse(str_detect(gene_emap, "sfp"), "Secondary metabolites", gene_label))
  # slice_sample(prop = 0.5, by = gene_label)

# Plot
# Pseudomonas log2FC across two contrasts (co vs mono (wt) and co vs mono (mt))
ps_p <- comps_wide_ps %>% 
  ggplot(aes(x = ps_covsmono_wt_l2fc, 
             y = ps_covsmono_mt_l2fc, 
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_ps, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_ps, !is.na(gene_label)),
             size = 3,
             pch = 21,
             col = "grey20") +
  coord_cartesian(
    xlim = c(-9, 9), 
    ylim = c(-4, 4)) +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  labs(x = "PS92 co (3610) - PS92 mono",
       y = "PS92 co (dhbA) - PS92 mono") +
  theme_bw() +
  theme(text = element_text(size = text_size),
        legend.position = "none"
  )


# Zoom of the previous plot
ps_zoom_p <- comps_wide_ps %>% 
  ggplot(aes(x = ps_covsmono_wt_l2fc, 
             y = ps_covsmono_mt_l2fc, 
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_ps, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_ps, !is.na(gene_label)),
             size = 3,
             pch = 21,
             col = "grey20") +
  coord_cartesian(
    xlim = c(-8, -lfc_thresh), 
    ylim = c(-lfc_thresh, lfc_thresh)) +
  geom_label_repel(data = dplyr::filter(comps_wide_ps, !is.na(gene_label),
                                        ps_covsmono_wt_l2fc < -lfc_thresh),
                   aes(label = gene),
                   label.size = 0.01,
                   max.overlaps = 7,
                   box.padding = 0.5,
                   nudge_y = 0.1,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20,
                   show.legend = F) +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  labs(x = "PS92 co (3610) - PS92 mono",
       y = "PS92 co (dhbA) - PS92 mono") +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    legend.position = "bottom"
  ) 

# Bacillus log2FC across two contrasts (co vs mono (wt) and co vs mono (mt))
covsmono_p <- comps_wide_bs %>% 
  ggplot(aes(x = bs_covsmono_l2fc, 
             y = mt_covsmono_l2fc,
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_bs, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_bs, !is.na(gene_label),
                                  gene_label != "Biofilm"),
             size = 3,
             pch = 21,
             col = "grey20") +
  labs(y = "dhbA co - dhbA mono (log2FC)",
       x = "3610 co - 3610 mono (log2FC)") +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "bottom"
  )

# Bacillus log2FC across two contrasts (mt vs wt (mono) and mt vs wt (co))
mtvswt_p <- comps_wide_bs %>% 
  ggplot(aes(x = bs_mtvswt_l2fc, 
             y = bs_co_mtvswt_l2fc,
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_bs, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_bs, !is.na(gene_label),
                                  gene_label != "Biofilm"),
             size = 3,
             pch = 21,
             col = "grey20") +
  labs(y = "dhbA co - 3610 co (log2FC)",
       x = "dhbA mono - 3610 mono (log2FC)") +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = text_size),
    text = element_text(size = text_size),
    legend.position = "none"
  )

# Combine
ps_zoom_p | ps_p | covsmono_p | mtvswt_p

ggsave("fig4_plots.pdf",
       path = here("figures", "figure_4_rnamsi"),
       height = 60,
       width = 185,
       units = "mm",
       scale = 2)


# Presentation
# Bacillus log2FC across two contrasts (co vs mono (wt) and co vs mono (mt))
covsmono_p <- comps_wide_bs %>% 
  ggplot(aes(x = bs_covsmono_l2fc, 
             y = mt_covsmono_l2fc,
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_bs, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_bs, !is.na(gene_label),
                                  gene_label != "Biofilm"),
             size = 3,
             pch = 21,
             col = "grey20") +
  labs(y = "dhbA co - dhbA mono (log2FC)",
       x = "3610 co - 3610 mono (log2FC)") +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 20),
    text = element_text(size = 20),
    legend.position = "bottom"
  )

# Bacillus log2FC across two contrasts (mt vs wt (mono) and mt vs wt (co))
mtvswt_p <- comps_wide_bs %>% 
  ggplot(aes(x = bs_mtvswt_l2fc, 
             y = bs_co_mtvswt_l2fc,
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_bs, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_bs, !is.na(gene_label),
                                  gene_label != "Biofilm"),
             size = 3,
             pch = 21,
             col = "grey20") +
  labs(y = "dhbA co - 3610 co (log2FC)",
       x = "dhbA mono - 3610 mono (log2FC)") +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 20),
    text = element_text(size = 20),
    legend.position = "bottom"
  )

plots_aligned <- align_patches(covsmono_p, mtvswt_p)

plots_aligned[[1]]
ggsave("fig4c.png",
       path = here("figures", "presentation"),
       dpi = 600,
       height = 100,
       width = 120,
       units = "mm",
       scale = 2)

plots_aligned[[2]]
ggsave("fig4d.png",
       path = here("figures", "presentation"),
       dpi = 600,
       height = 100,
       width = 120,
       units = "mm",
       scale = 2)


ps_p <- comps_wide_ps %>% 
  ggplot(aes(x = ps_covsmono_wt_l2fc, 
             y = ps_covsmono_mt_l2fc, 
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_ps, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_ps, !is.na(gene_label)),
             size = 3,
             pch = 21,
             col = "grey20") +
  coord_cartesian(
    xlim = c(-9, 9), 
    ylim = c(-4, 4)) +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  labs(x = "PS92 co (3610) - PS92 mono",
       y = "PS92 co (dhbA) - PS92 mono") +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "bottom"
  )


# Zoom of the previous plot
ps_zoom_p <- comps_wide_ps %>% 
  ggplot(aes(x = ps_covsmono_wt_l2fc, 
             y = ps_covsmono_mt_l2fc, 
             fill = gene_label)) +
  geom_hline(yintercept = lfc_thresh, linetype = 2) +
  geom_hline(yintercept = -lfc_thresh, linetype = 2) +
  geom_vline(xintercept = lfc_thresh, linetype = 2) +
  geom_vline(xintercept = -lfc_thresh, linetype = 2) +
  geom_point(data = dplyr::filter(comps_wide_ps, is.na(gene_label)),
             size = 0.5,
             alpha = 0.5,
             col = "grey50") +
  geom_point(data = dplyr::filter(comps_wide_ps, !is.na(gene_label)),
             size = 3,
             pch = 21,
             col = "grey20") +
  coord_cartesian(
    xlim = c(-8, -lfc_thresh), 
    ylim = c(-lfc_thresh, lfc_thresh)) +
  geom_label_repel(data = dplyr::filter(comps_wide_ps, !is.na(gene_label),
                                        ps_covsmono_wt_l2fc < -lfc_thresh),
                   aes(label = gene),
                   label.size = 0.01,
                   max.overlaps = 7,
                   box.padding = 0.5,
                   nudge_y = 0.1,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20,
                   show.legend = F) +
  scale_fill_brewer(palette = "Set2",
                    name = NULL) +
  labs(x = "PS92 co (3610) - PS92 mono",
       y = "PS92 co (dhbA) - PS92 mono") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = "bottom"
  ) 

plots_aligned <- align_patches(ps_p, ps_zoom_p)

plots_aligned[[1]]
ggsave("fig4b.png",
       path = here("figures", "presentation"),
       dpi = 600,
       height = 100,
       width = 120,
       units = "mm",
       scale = 2)

plots_aligned[[2]]
ggsave("fig4a.png",
       path = here("figures", "presentation"),
       dpi = 600,
       height = 100,
       width = 120,
       units = "mm",
       scale = 2)



#### Figure S4: RNA-seq QC plots
#### Figure S4: Transcriptome QC and functional analysis
