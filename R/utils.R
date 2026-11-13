

ss_features <- c(setNames("Alpha helix (4-12)", "H"),
                 setNames("Isolated beta-bridge residue", "B"),
                 setNames("Strand", "E"),
                 setNames("3-10 helix", "G"),
                 setNames("Pi helix", "I"),
                 setNames("Turn", "T"),
                 setNames("Bend", "S"),
                 setNames("Kappa helix", "P"))


uniprotID_to_uniprotName <- function(x) {

  id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

  converter <- setNames(id_map[["Entry Name"]], id_map[["Entry"]])

  unname(converter[x])

}

tibble_split <- function(..., names = NULL) {
  tmp <- stringr::str_split(..., simplify = TRUE)
  if(is.null(names)) {
    names <- paste0("split", 1:ncol(tmp))
  }
  if(length(names) != ncol(tmp)) {
    stop("names not equal to number of split fields")
  }
  colnames(tmp) <- names
  tmp %>% as_tibble
}

tibble_table <- function(x) {
  tmp <- table(x)
  tibble(name = names(tmp),
         freq = as.vector(unname(tmp)))
}


extra_color <- c("#FD7446FF","#FD8CC1FF","#197EC0FF","#FED439FF",
                 "#D5E4A2FF", "#D2AF81FF", "#C80813FF", "#46732EFF", "#71D0F5FF",
                 "#370335FF", "#075149FF", "#91331FFF", "#1A9993FF", "#709AE1FF", "#FD744680", "#46732E80", "#FF00FF80", "#197EC080")

AA_colors <- c(A = "#FD744680",
               I = "#F05C3BAA",
               L = "#C8081380",
               M = "#FF6666",
               `F` = "#FF6961",
               W = "#FA8072",
               V = "#FFB6C1",
               S = "#1A999380",
               `T` = "#1A9993AA",
               N = "#46732E80",
               Q = "#D5E4A2AA",
               K = "#709AE180",
               R = "#71D0F5FF",
               H = "#197EC080",
               D = "#C7158580",
               E = "#FF00FF80",
               G = "#FFA500AA",
               P = "#D2AF81CC",
               Y = "#00FFFF80",
               C = "#FED439CC",
               `-` = "#8A919760")

factor_to_uniprotFeature <- function(tf,
                                     source = "label1",
                                     evidence = "label2") {

  uni_tf <- unique(tf)

  bind_rows(
    lapply(uni_tf[!uni_tf == "-"], \(x) {
      vec <- which(tf == x)

      breaks <- c(0, which(diff(vec) != 1), length(vec))

      sequences <- lapply(seq_along(breaks[-1]), \(i) vec[(breaks[i] + 1):breaks[i + 1]])

      bind_rows(lapply(sequences, \(y) tibble(type = x,
                                              evidence = evidence,
                                              start = min(y),
                                              end = max(y),
                                              source = source)
      ))

    })
  )

}

distance <- function(s, p) {
  sqrt(sum((s - p)^2))
}

get_mode <- function(x) {
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}


yo <- function(text = "something happened with some code") {
  temp.file <- paste0(tempdir(), "/error.txt")
  file.create(temp.file)
  writeLines(text, con = temp.file)
  system(paste("sendmail", "kbrulois@stanford.edu", "<", temp.file))
  system(paste("sendmail", "5628226527@vtext.com", "<", temp.file))
}
