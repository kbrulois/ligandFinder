

# after new6 finishes scoring:
saveRDS(nn_input_comb, "~/AF2_analysis/nn_input_comb.rds")
saveRDS(nn_input,      "~/AF2_analysis/nn_input.rds")        # optional: train/val membership, val diagnostics

# the trained models are keras/python objects -> NOT saveRDS. Only needed if you
# want to re-SCORE new windows later (not for plotting):
options(error = NULL)   # in case recover is still armed

for (tm in names(models)) {
  keras3::save_model(models[[tm]], path.expand(sprintf("~/AF2_analysis/model_%s.keras", tm)))
}

saveRDS(mget(c("class_cols","classes","real_cols","seq_len","n_channels","all_params3")),
        "~/AF2_analysis/plot_ctx.rds")
