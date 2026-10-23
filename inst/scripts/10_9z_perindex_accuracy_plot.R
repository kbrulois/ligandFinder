suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
meta <- tribble(~tag, ~label, ~anchor,
  "C_near", "C terminus", 28, "C_ctpad", "C terminus, CT-padded", 28,
  "N_near", "N terminus", 8,  "N_ntpad", "N terminus, NT-padded", 8)
cols <- c("C terminus" = "#D95F0E", "C terminus, CT-padded" = "#8C2D04",
          "N terminus" = "#2BA84A", "N terminus, NT-padded" = "#006D2C")
nice <- c(CT_cleavage_context = "CT context", NT_cleavage_context = "NT context",
          pep_other = "pep_other", pep_pocket = "pep_pocket", none = "none", padding = "padding")
J <- map(set_names(meta$tag, meta$label), ~ fromJSON(file.path(S2, paste0("pix_", .x, ".json")),
                                                     simplifyVector = FALSE))
pos <- imap_dfr(J, function(j, lab) map_dfr(j$by_position, function(p)
  tibble(label = lab, pos = p$pos, acc = p$acc %||% NA_real_, n = p$n)))
anc <- meta %>% select(label, anchor)
p1 <- ggplot(pos, aes(pos, acc, colour = label)) +
  geom_vline(data = anc, aes(xintercept = anchor), linetype = 2, colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 0.8) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = c(1, 8, 15, 22, 28, 36)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Per-residue accuracy along the window",
       subtitle = "residues whose true label is not padding; dashed = the anchor (28 for C, 8 for N)",
       x = "position in window", y = "correct class")

pc <- imap_dfr(J, function(j, lab) imap_dfr(j$per_class, function(d, nm)
  if (is.null(d$recall)) NULL else tibble(label = lab, class = nice[[nm]],
                                          recall = d$recall, n = d$n)))
pc <- pc %>% mutate(class = factor(class, levels = rev(c("CT context", "NT context",
                                                         "pep_other", "pep_pocket", "none", "padding"))))
p2 <- ggplot(pc, aes(recall, class, colour = label)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Recall by class", subtitle = "each trunk is best at its own terminus's context",
       x = "recall on held-out windows", y = NULL)

cm <- imap_dfr(J, function(j, lab) {
  m <- do.call(rbind, lapply(j$confusion, unlist)); cls <- unlist(j$classes)
  dimnames(m) <- list(true = cls, pred = cls)
  as_tibble(as.table(m), n = "n") %>%
    group_by(true) %>% mutate(frac = n / sum(n)) %>% ungroup() %>% mutate(label = lab)
})
lv <- c("CT_cleavage_context", "NT_cleavage_context", "pep_other", "pep_pocket", "none", "padding")
cm <- cm %>% mutate(true = factor(nice[as.character(true)], levels = rev(nice[lv])),
                    pred = factor(nice[as.character(pred)], levels = nice[lv]))
p3 <- ggplot(cm, aes(pred, true, fill = frac)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(frac >= 0.08, sprintf("%.0f", 100 * frac), "")),
            size = 2.5, colour = "grey10") +
  facet_wrap(~ label, nrow = 1) +
  scale_fill_gradient(low = "#F4F7F9", high = "#197EC0", labels = scales::percent, name = "of true class") +
  labs(title = "Where the errors go", subtitle = "row-normalised: each row is one true class, % predicted as each class",
       x = "predicted", y = "true") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.4, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA))
hdr <- map_chr(names(J), ~ sprintf("%s: %.0f%% overall, %.0f%% on the four real classes",
                                   .x, 100 * J[[.x]]$acc_all, 100 * J[[.x]]$acc_real_only))
p <- ((p1 & th) | (p2 & th)) / (p3 & th + theme(axis.text.x = element_text(angle = 40, hjust = 1))) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Per-index head: how well each residue of a held-out window is classified",
    subtitle = paste0("`none` is 70% of residues (every negative window is all `none`), so overall accuracy is mostly that class.\n",
                      paste(hdr, collapse = "   ·   ")),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7.6, colour = "grey30")))
f <- file.path(A2, "pepend_perindex_accuracy.svg")
ggsave(f, p, width = 13.5, height = 8.2)
write_csv(pc, file.path(A2, "pepend_perindex_recall.csv"))
message("wrote ", f)
