out_dir <- "/dskh/nobackup/cliang/TRACE/fig4"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(patchwork)
})

poag_ids <- c("D011", "D023")
p_threshold <- 0.01
dataset_labels <- c(
  D011 = "D011 POAG blood",
  D023 = "D023 POAG organoid"
)

files <- file.path(out_dir, paste0(poag_ids, "_asgard_drug_results_long.csv"))
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing ASGARD result files:\n", paste(missing_files, collapse = "\n"))
}

drug_summary <- bind_rows(lapply(files, read_csv, show_col_types = FALSE)) |>
  filter(dataset_id %in% poag_ids) |>
  group_by(dataset_id, dataset_name, Drug.name) |>
  summarise(
    min_P = min(P.value, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    best_rank = min(rank_within_cluster, na.rm = TRUE),
    n_passing_clusters = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    support = case_when(
      min_FDR < 0.05 ~ "FDR significant",
      min_P < p_threshold ~ "Nominal support",
      TRUE ~ "No support"
    ),
    supported = support != "No support"
  )

selected_drug_table <- drug_summary |>
  filter(supported) |>
  group_by(Drug.name) |>
  summarise(
    n_supported = n_distinct(dataset_id),
    n_significant = n_distinct(dataset_id[min_FDR < 0.05]),
    best_any_fdr = min(min_FDR, na.rm = TRUE),
    best_any_p = min(min_P, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_supported), desc(n_significant), best_any_fdr, best_any_p) |>
  mutate(include_priority = n_significant > 0)

selected_drugs <- bind_rows(
  selected_drug_table |> filter(include_priority),
  selected_drug_table |> filter(!include_priority) |> slice_head(n = 28)
) |>
  arrange(desc(n_significant), desc(n_supported), best_any_fdr, best_any_p)

if (nrow(selected_drugs) == 0) {
  stop("No supported POAG drugs found.")
}

lead_drugs <- c("belinostat", "tranylcypromine")
selected_drugs <- selected_drugs |>
  mutate(lead_order = match(Drug.name, lead_drugs)) |>
  arrange(if_else(is.na(lead_order), 3L, lead_order), desc(n_significant),
          desc(n_supported), best_any_fdr, best_any_p) |>
  select(-lead_order)

drug_labels <- setNames(tools::toTitleCase(selected_drugs$Drug.name), selected_drugs$Drug.name)
drug_labels[["belinostat"]] <- "Belinostat (HDAC inhibitor)"
drug_labels[["tranylcypromine"]] <- "Tranylcypromine (LSD1 inhibitor)"

heatmap_data <- drug_summary |>
  filter(Drug.name %in% selected_drugs$Drug.name) |>
  mutate(
    dataset_label = factor(dataset_labels[dataset_id], levels = dataset_labels[poag_ids]),
    Drug.name = factor(Drug.name, levels = selected_drugs$Drug.name),
    # Rescale the visually relevant leading ranks only: rank 1 = 1 and
    # rank 10 = 0. Ranks below the displayed top-ten range remain at 0.
    rank_score = pmax(0, 1 - ((best_rank - 1) / 9)),
    is_lead = best_rank == 1
  )

heatmap <- ggplot(heatmap_data, aes(x = Drug.name, y = dataset_label, fill = rank_score)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(
    data = \(d) dplyr::filter(d, is_lead),
    label = "*", color = "white", size = 8, fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "#fde0dd", high = "#a50f15", limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = "Normalized\nrank score"
  ) +
  scale_x_discrete(labels = drug_labels) +
  labs(
    title = "Dataset-specific POAG ASGARD drug rankings",
    subtitle = "Colour bar: rank 1 = 1; ranks 2-10 are scaled to 0-1; ranks >10 = 0",
    caption = "* Highest-ranked ASGARD candidate within its dataset",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.caption = element_text(hjust = 0, face = "italic", size = 9)
  )

drug_membership <- drug_summary |>
  filter(supported) |>
  distinct(Drug.name, dataset_id) |>
  group_by(Drug.name) |>
  summarise(
    category = case_when(
      all(poag_ids %in% dataset_id) ~ "Shared",
      "D011" %in% dataset_id ~ "D011 only",
      "D023" %in% dataset_id ~ "D023 only",
      TRUE ~ "Other"
    ),
    .groups = "drop"
  )

overlap_summary <- drug_membership |>
  count(category, name = "n_drugs") |>
  mutate(
    category = factor(category, levels = c("Shared", "D011 only", "D023 only"))
  )

barplot <- ggplot(overlap_summary, aes(x = n_drugs, y = category, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_drugs), hjust = -0.25, size = 3.5) +
  scale_fill_manual(
    values = c("Shared" = "#8c2d04", "D011 only" = "#2b6cb0", "D023 only" = "#de2d26"),
    guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Supported Drug Recurrence",
    subtitle = "P < 0.01",
    x = "Number of drugs",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.ticks.y = element_blank()
  )

combined <- heatmap / barplot + plot_layout(heights = c(1.15, 0.85))

ggsave(file.path(out_dir, "poag_asgard_drug_support_and_recurrence.pdf"), combined, width = 9.6, height = 8.2)
ggsave(file.path(out_dir, "poag_asgard_drug_support_and_recurrence.png"), combined, width = 9.6, height = 8.2, dpi = 300)

# Alternative view: an ordered lollipop plot makes exact rank differences readable.
lollipop_data <- drug_summary |>
  group_by(dataset_id) |>
  arrange(best_rank, min_FDR, Drug.name, .by_group = TRUE) |>
  slice_head(n = 10) |>
  mutate(plot_order = row_number()) |>
  ungroup() |>
  mutate(
    # After repeated D011 Belinostat hits are collapsed, use consecutive
    # positions for its distinct displayed drugs.  Keep D023's original
    # within-cluster positions unchanged.
    display_rank = case_when(
      dataset_id == "D011" & Drug.name %in%
        c("alverine", "benzonatate", "chloroxine", "dobutamine") ~ 7L,
      dataset_id == "D011" ~ plot_order,
      TRUE ~ best_rank
    ),
    dataset_label = factor(dataset_labels[dataset_id], levels = dataset_labels[poag_ids]),
    display_drug = case_when(
      Drug.name == "belinostat" ~ "Belinostat (HDAC inhibitor)",
      Drug.name == "tranylcypromine" ~ "Tranylcypromine (LSD1 inhibitor)",
      TRUE ~ tools::toTitleCase(Drug.name)
    )
  )

lollipop_plot <- ggplot(lollipop_data, aes(x = display_rank, y = -plot_order, colour = dataset_label)) +
  geom_vline(xintercept = 1:10, colour = "grey90", linewidth = 0.35) +
  geom_segment(aes(x = 10, xend = display_rank, yend = -plot_order),
               colour = "grey70", linewidth = 0.75) +
  geom_point(size = 5.2) +
  geom_text(aes(x = 13.8, label = display_drug, fontface = if_else(display_rank == 1, "bold", "plain")),
            colour = "black", hjust = 0, size = 3.6) +
  facet_wrap(~dataset_label, ncol = 1) +
  scale_colour_manual(
    values = c("D011 POAG blood" = "#2b6cb0", "D023 POAG organoid" = "#c61c1c"),
    guide = "none"
  ) +
  scale_x_reverse(
    breaks = 10:1, limits = c(14.2, 0.4),
    name = "ASGARD within-cluster rank (1 = highest)"
  ) +
  scale_y_continuous(breaks = NULL) +
  labs(
    title = "POAG ASGARD drug rankings",
    subtitle = "Top 10 candidates per dataset; rank 1 is the leading candidate",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.ticks.y = element_blank(),
    panel.spacing = grid::unit(0.8, "lines")
  )

ggsave(file.path(out_dir, "poag_asgard_drug_ranking_lollipop.pdf"), lollipop_plot, width = 9.6, height = 7.6)
ggsave(file.path(out_dir, "poag_asgard_drug_ranking_lollipop.png"), lollipop_plot, width = 9.6, height = 7.6, dpi = 300)

# Main Figure 4: show the two dataset-specific lollipop panels alongside the
# recurrence summary in a single row.
lollipop_d011 <- (lollipop_plot %+% filter(lollipop_data, dataset_id == "D011")) +
  labs(title = NULL, subtitle = NULL, x = NULL) +
  theme(plot.margin = margin(5.5, 2, 5.5, 2))
lollipop_d023 <- (lollipop_plot %+% filter(lollipop_data, dataset_id == "D023")) +
  labs(title = NULL, subtitle = NULL, x = NULL) +
  theme(plot.margin = margin(5.5, 2, 5.5, 2))

three_column_figure <- lollipop_d011 + lollipop_d023 + barplot +
  plot_layout(widths = c(1.25, 1.25, 1)) +
  plot_annotation(
    title = "POAG ASGARD drug rankings and recurrence",
    subtitle = "Lollipop panels show the top 10 candidates per dataset; rank 1 is highest"
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

ggsave(file.path(out_dir, "poag_asgard_drug_support_and_recurrence.pdf"),
       three_column_figure, width = 14.5, height = 5.6)
ggsave(file.path(out_dir, "poag_asgard_drug_support_and_recurrence.png"),
       three_column_figure, width = 14.5, height = 5.6, dpi = 300)

write_csv(heatmap_data, file.path(out_dir, "poag_asgard_drug_support_heatmap_data.csv"))
write_csv(lollipop_data, file.path(out_dir, "poag_asgard_drug_ranking_lollipop_data.csv"))
write_csv(overlap_summary, file.path(out_dir, "poag_asgard_drug_recurrence_summary.csv"))

message("Wrote POAG ASGARD support and recurrence figure to ", out_dir)
