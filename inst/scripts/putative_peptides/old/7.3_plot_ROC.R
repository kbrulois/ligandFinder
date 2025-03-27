

extra_color <- c("#FD7446FF","#FD8CC1FF","#FED439FF", "#197EC0FF", "#46732EFF", "#C80813FF",
                 "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF",
                 "#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF")



library(ggplot2)
library(pROC)


set.seed(1)
secretome$known_ligand2 <- ifelse(secretome$dbr_type == "gtp", 1, 0)
samp <- sample(c(TRUE, FALSE), nrow(secretome), replace=TRUE, prob=c(0.7,0.3))
train <- secretome[samp, ]
test <- secretome[!samp, ]

secretome_aa <- secretome_aa %>%
  mutate(conservation = scales::rescale(log(conservation + 0.1), c(1,0)))

samp <- sample(c(TRUE, FALSE), length(unique(secretome_aa$accession)), replace=TRUE, prob=c(0.7,0.3))
train <- secretome_aa[secretome_aa[["accession"]] %in% unique(secretome_aa[["accession"]])[samp], ]
test <- secretome_aa[secretome_aa[["accession"]] %in% unique(secretome_aa[["accession"]])[!samp], ]



models <- list(conservation = "conservation",
               relASA = "relASA",
               pathogenicity = "pathogenicity",
               `conservation +\nrelASA +\npathogenicity` = c("conservation", "relASA", "pathogenicity"))

models2 <- list(conservation = "prox_cons_3",
                relASA = "prox_relASA_3",
                pathogenicity = "prox_afm_3",
                score = "score_prox",
                `conservation +\nrelASA +\npathogenicity` = c("prox_cons_3", "prox_relASA_3", "prox_afm_3"))

model_sets <- list(entire = models,
                   proximal = models2)

plots <- lapply(names(model_sets)[1], \(y) {
  
  
  modelz <- lapply(names(model_sets[[y]]), \(x) {
    message("computing ", y, x)
    
    formula <- as.formula(paste("known_peptide", "~", paste(model_sets[[y]][[x]], collapse = " + ")))
    
    
    model_glm <- glm(formula, family="binomial", data=train)
    
    predicted_glm <- predict(model_glm, test, type="response")
    
    glm = roc(test$known_peptide, predicted_glm, ci = TRUE)
    
    
    if(x == "conservation +\nrelASA +\npathogenicity") {
      
      train2 <- train %>% select(known_peptide, conservation, relASA, pathogenicity) %>% drop_na()
      train3 <- train2 %>% select(-known_peptide, conservation, relASA, pathogenicity)
      test2 <- test %>% select(known_peptide, conservation, relASA, pathogenicity) %>% drop_na()
      
      message("computing neural net")
      
      model_nn = neuralnet::neuralnet(formula,
                                      data=train2,
                                      hidden=3,
                                      linear.output = FALSE
      )
      
      message("computing random forest")
      
      
      # model_rf <- randomForest::randomForest(x = train3, 
      #                                        y = train2$known_peptide)
      
      predicted_nn <- predict(model_nn, test, type = "response")
      #predicted_rf <- predict(model_rf, test2 %>% select(-known_peptide), type = "response")
      
      nn = roc(test$known_peptide, predicted_nn, ci = TRUE)
      #rf = roc(test2$known_peptide, predicted_rf, ci = TRUE)
      
      to_return <- list(glm = glm, nn = nn)
    } else {
      to_return <- list(glm = glm)
    }
    
    names(to_return) <- paste0(x, "_", names(to_return))
    return(to_return)
  })
  
  modelz <- do.call(c, modelz)
  
  modelz <- modelz[!sapply(modelz, is.null)]
  
  names(modelz) <- paste(names(modelz), "\n[AUC:", round(sapply(modelz, auc), 4), "]")
  
  
  ggroc(modelz, linewidth = 1) + ggplot2::theme_bw() + scale_color_discrete(name = "", type = extra_color) +
    ggtitle(paste("ROC Analysis"), subtitle = y) +  theme(
      legend.key.height = unit(3, "lines")  
    )
  
})


p_al <- Reduce(`+`, plots) +
  
  patchwork::plot_layout(ncol = 2)


ggsave("~/peptide_alg/ROC_analysis_per_aa_w_conservation_transformed.svg", plot = p_al, device = svglite::svglite, width = 12, height = 5)



