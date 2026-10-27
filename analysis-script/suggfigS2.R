# Reads results/uncertainty_per_cell.csv, written by fig4_uncertainty.R; run that script first.
suppressMessages({library(ggplot2); library(ggrepel); library(dplyr)})
out_dir <- "results"
# Cell type of each published cluster, from Supplementary Table 1 of Youkilis & Bassnett (2021)
s1 <- data.frame(cluster = 1:22, celltype = c("Rods", "U - Fibroblast", "Amacrine", "Macrophage", "Cone Bipolar", "Melanocyte", "M\u00fcller Glia", "RPE",
  "PCE", "S - Fibroblast", "Rod Bipolar", "NPCE", "Cones", "Vascular Endothelium", "rM Macrophage", "N/K T-Cell", "Erythrocyte", "Smooth Muscle",
  "Neutrophil", "Schwann Cell", "B - Cell", "RGC"))
types <- s1$celltype
x <- read.csv(file.path(out_dir, "uncertainty_per_cell.csv"), stringsAsFactors = FALSE, encoding = "UTF-8"); x$type <- factor(x$type, levels = types)
stopifnot(nrow(x) == 10024, !anyNA(x$type))

ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"
type_cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948", "#86b6ef", "#f4a582", "#7fd1b3", "#f6cf6a",
                        "#b05c84", "#66b266", "#9085e9", "#a83232", "#184f95", "#a84a1f", "#0f7a55", "#8a6400", "#6d6c66", "#b9b8b0"), types)
bs <- 7; txt <- (bs - 1) / ggplot2::.pt
th <- theme_minimal(base_size = bs) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(colour = ink2),
  axis.text = element_text(colour = ink, size = bs - 0.5), axis.title = element_text(colour = ink2, size = bs),
  plot.title = element_text(colour = ink, face = "bold", size = bs + 1, margin = margin(b = 2)), plot.title.position = "plot",
  plot.subtitle = element_text(colour = ink2, size = bs - 0.3, lineheight = 1.1, margin = margin(b = 4)), legend.text = element_text(size = bs - 0.3),
  legend.title = element_text(size = bs - 0.3), plot.margin = margin(5, 7, 5, 5), plot.background = element_rect(fill = "white", colour = NA))
box <- theme(panel.border = element_rect(colour = ink2, fill = NA, linewidth = 0.4))

# Supplementary figure: cell type proportion against uncertainty
d <- x %>% group_by(type) %>% summarise(n = n(), uncertainty = mean(uncertainty), .groups = "drop") %>% mutate(pct = 100 * n / sum(n))
ct <- cor.test(d$pct, d$uncertainty, method = "spearman", exact = FALSE)        # rank-based: proportions are very skewed
cat(sprintf("\ncell type proportion vs mean uncertainty: Spearman rho = %.2f, p = %.2f\n", ct$estimate, ct$p.value))
supp <- ggplot(d, aes(pct, uncertainty)) + geom_point(aes(fill = type), shape = 21, size = 2.4, colour = "white", stroke = 0.5) +
  geom_text_repel(aes(label = type), colour = ink, size = txt * 0.95, segment.colour = muted, segment.size = 0.2, min.segment.length = 0.2, max.overlaps = Inf, seed = 3, box.padding = 0.25) +
  scale_fill_manual(values = type_cols, guide = "none") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 25, 5), expand = expansion(mult = c(0.02, 0.06)), labels = function(v) paste0(v, "%")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)), labels = function(v) paste0(v, "%")) +
  labs(title = "Cell type proportion against uncertainty", x = "Cell type proportion (% of all cells)", y = "Mean uncertainty") + th + box + theme(aspect.ratio = 1)

ggsave(file.path(out_dir, "supp_fig_proportion_vs_uncertainty.pdf"), supp, width = 3.8, height = 4.1, useDingbats = FALSE)
cat("\nwrote supp_fig_proportion_vs_uncertainty.pdf in", out_dir, "\n")
