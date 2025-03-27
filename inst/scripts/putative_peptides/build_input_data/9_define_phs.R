
secretome <- readRDS(paste0(s_localDir, "/processed/secretome_5.rds"))

all_features <- bind_rows(lapply(secretome$features, 
                                 \(x) {tibble(type = x[["type"]],
                                              source = x[["source"]])
                                 })) %>%
  distinct(type, .keep_all = TRUE)


sum(sapply(secretome$features, \(x) "alpha fold" %in% x[["source"]]))




secretome <- secretome %>%
  mutate(features = map(features, \(x) {
    x %>%
      mutate(start = as.integer(start),
             end = as.integer(end))
  })) %>%
  mutate(features = pmap(list(features, sequence_uni),
                         .f = \(features, sequence_uni) {
                           
                           Dibasic <- stringr::str_locate_all(sequence_uni, "KK|KR|RR|RK")[[1]]
                           C_sites <- stringr::str_locate_all(sequence_uni, "C")[[1]]
                           WY_sites <- stringr::str_locate_all(sequence_uni, "W|Y")[[1]]
                           num_DB <- nrow(Dibasic)
                           num_C <- nrow(C_sites)
                           num_WY <- nrow(WY_sites)
                           if(num_DB > 0) {
                             features <- bind_rows(features,
                                                   tibble(type = rep("Dibasic", num_DB),
                                                          evidence = rep("pattern search", num_DB),
                                                          start = Dibasic[,1],
                                                          end = Dibasic[,2],
                                                          source = "sites")
                             )
                           }
                           if(num_C > 0) {
                             features <- bind_rows(features,
                                                   tibble(type = rep("Cysteine", num_C),
                                                          evidence = rep("pattern search", num_C),
                                                          start = C_sites[,1],
                                                          end = C_sites[,2],
                                                          source = "sites")
                             )
                           }
                           if(num_WY > 0) {
                             features <- bind_rows(features,
                                                   tibble(type = rep("W or Y", num_WY),
                                                          evidence = rep("pattern search", num_WY),
                                                          start = WY_sites[,1],
                                                          end = WY_sites[,2],
                                                          source = "sites")
                             )
                           }
                           return(features)
                         })) %>%
  mutate(features = map(features,
                        .f = \(features) {
                          if(!"signal peptide" %in% features[["type"]]) {
                            
                            return(bind_rows(features,
                                             tibble(type = "signal peptide",
                                                    evidence = "F21",
                                                    start = 1,
                                                    end = 21,
                                                    source = "manual")
                            ))
                            
                          } else if(features %>% 
                                    filter(type == "signal peptide") %>% 
                                    pull(end) %>% 
                                    is.na) {
                            features[["end"]][features[["type"]] == "signal peptide"] <- 21
                            features[["source"]][features[["type"]] == "signal peptide"] <- "manual"
                            
                            return(features)
                          } else {
                            return(features)
                          }
                          
                        })) %>%
  mutate(features = pmap(list(features, sequence_uni),
                         .f = \(features, sequence_uni) {
                           ct <- nchar(sequence_uni)
                           bind_rows(features,
                                     tibble(type = "c-terminus",
                                            evidence = "lastAA",
                                            start = ct,
                                            end = ct,
                                            source = "manual")
                           )
                         })) %>%
  mutate(features = map(features, \(x) x %>%
                          drop_na(type)))






###define putative peptide ligands (ppl)

summary(sapply(secretome$features, \(x) sum(x[["type"]] == "Dibasic")))  

x <- tibble(type = c("signal peptide",
                     "helix",
                     "Dibasic", #very close to start
                     "Dibasic", #truncated by sp and beta
                     "Dibasic", #KK_KR in offlimits region
                     "strand",
                     "Dibasic",#free
                     "Dibasic", #in range of next KK_KR
                     "Dibasic",
                     "c-terminus"),#truncated by c terminus
            evidence = 1,
            start = c(1,24,22,36,52,51,85,120,130,145),
            end = c(21,31,23,37,53,58,86,121,131,145)
)



define_dbr <- function(x, 
                       offlimits = c("signal peptide", "signal peptide (manual)", "E"),
                       offlimits_alt = c("signal peptide", "signal peptide (manual)", "strand"),
                       search_width = 10) {
  
  if(!"Dibasic" %in% x[["type"]]) {
    return(x)
  } else {
    if(!"alpha fold" %in% x[["source"]]) {
      offlimits <- offlimits_alt
    }
    kk_kr <- which(x[["type"]] == "Dibasic")
    ol_features <- x %>% filter(type %in% offlimits)
    
    bind_rows(
      lapply(kk_kr, \(y) {
        site <- unlist(x[y,c("start", "end")], use.names = FALSE)
        
        cond <- ol_features %>%
          mutate(overlap = start <= site[2] & end >= site[1]) %>%
          summarize(overlap = any(overlap)) %>%
          pull(overlap)
        
        if(cond) {
          return(tibble(type = character(), 
                        evidence = character(), 
                        start = integer(), 
                        end = integer(),
                        source = character()))
        } else {
          
          roi <- tibble(type = c(paste("phs_dbr", paste0(site[1], "-", site[2]), "us", sep = "_"),
                                 paste("phs_dbr", paste0(site[1], "-", site[2]), "ds", sep = "_")),
                        evidence = c("phs_dbr", "phs_dbr"),
                        start = c(site[1] - (search_width - 1), site[2]),
                        end = c(site[1], site[2] + search_width - 1),
                        source = "phs_dbr")
          
          
          #upstream
          site <- unlist(roi[1,c("start", "end")], use.names = FALSE)
          
          if(site[1] <= 0) site[1] <- roi[["start"]][1] <- 1
          
          trim <- x %>%
            filter(type %in% c(offlimits, "Dibasic") & end >= site[1] & end <= site[2]) 
          
          if(nrow(trim) > 0) {
            roi[["start"]][1] <- trim %>%
              summarize(trim = max(end)) %>%
              pull(trim)
          }
          
          #downstream
          site <- unlist(roi[2,c("start", "end")], use.names = FALSE)
          
          cterm <- x[["end"]][x[["type"]] == "c-terminus"]
          if(site[2] > cterm) site[2] <- roi[["end"]][2] <- cterm
          
          trim <- x %>%
            filter(type %in% c(offlimits, "Dibasic") & start <= site[2] & start >= site[1])
          
          if(nrow(trim) > 0) {
            roi[["end"]][2] <- trim %>%
              summarize(trim = min(start)) %>%
              pull(trim)
          }         
          
          return(roi)
          
        }
        
      })
    ) -> rois
    
    rois <- rois %>%
      distinct(start, end, .keep_all = TRUE) %>%
      filter(!start >= end)
    
    
    return(bind_rows(rois, x))
  }
}

define_termini <- function(x, 
                           search_width = 10) {
  
  start_n <- x %>% 
    filter(type == "signal peptide") %>%
    mutate(end = end + 1) %>%
    pull(end)
  
  end_c <- x %>% 
    filter(type == "c-terminus") %>%
    pull(start)
  
  end_n_limit <- start_n + (end_c - start_n)/3
  
  start_c_limit <- start_n + 2 * (end_c - start_n)/3
  
  if(start_n < end_c) {
    rois <- tibble(type = c(paste("phs_c-term_", paste0(max(end_c - (search_width - 1), start_n), "-", end_c), "us", sep = "_"),
                            paste("phs_n-term_", paste0(start_n, "-", min(start_n + (search_width - 1), end_c)), "ds", sep = "_")),
                   evidence = c("phs_c-term", "phs_n-term"),
                   start = c(max(end_c - (search_width - 1), start_n), start_n),
                   end = c(end_c, min(start_n + (search_width - 1), end_c)),
                   source = c("phs_c-term", "phs_n-term"))
  } else {
    rois <- tibble(type = character(), 
                   evidence = character(), 
                   start = integer(), 
                   end = integer(), 
                   source = character(),
                   .rows = 0)
  }
  
  num_n_cys <- x %>%
    filter(type == "Cysteine" & between(start, start_n + search_width, max(start_n + search_width, end_n_limit)))
  
  if(nrow(num_n_cys) > 0) {
    
    end_n <- num_n_cys %>% 
      summarize(min_start = min(start)) %>% 
      pull(min_start)
    
    rois_cys_n <- tibble(type = c(paste("phs_n-term-Cys_", paste0(max(end_n - (search_width + 1), start_n), "-", end_n), "us", sep = "_")),
                         evidence = "n-term-Cys",
                         start = max(end_n - (search_width + 1), start_n),
                         end = end_n,
                         source = "phs_n-term-Cys")
  } else {
    rois_cys_n <- tibble(type = character(), 
                         evidence = character(), 
                         start = integer(), 
                         end = integer(), 
                         source = character(),
                         .rows = 0)
  }
  
  num_c_cys <- x %>%
    filter(type == "Cysteine" & between(start, min(start_c_limit , end_c - (search_width + 1)), end_c - (search_width + 1)))
  
  if(nrow(num_c_cys) > 0) {
    
    start_c <- num_c_cys %>% 
      summarize(max_start = max(start)) %>% 
      pull(max_start)
    start_c <- min(start_c, end_c - (search_width - 1))
    
    rois_cys_c <- tibble(type = c(paste("phs_c-term-Cys_", paste0(min(end_c - (search_width - 1), start_c), "-", end_c), "us", sep = "_")),
                         evidence = "c-term-Cys",
                         start = start_c,
                         end = min(start_c + (search_width - 1), end_c),
                         source = "phs_c-term-Cys")
  } else {
    rois_cys_c <- tibble(type = character(), 
                         evidence = character(), 
                         start = integer(), 
                         end = integer(), 
                         source = character(),
                         .rows = 0)
  }
  
  return(bind_rows(rois, rois_cys_n, rois_cys_c, x))
  
}

secretome <- secretome %>%
  mutate(features = map(features, .f = define_dbr)) %>%
  mutate(features = map(features, .f = define_termini))


add_site_stats <- function(features, cons, dssp, sequence_uni) {
  
  regions <- grep("^dbr_", features[["type"]], value = TRUE)
  
  regions2 <- features %>%
    filter(source == "gtp") %>%
    pull(type)
  
  regions <- c(regions, regions2)
  
  if(length(regions) == 0) {
    return(features)
  } else {
    
    aa_sequence <- strsplit(sequence_uni, "")[[1]]
    
    if("frequency" %in% names(cons$ms)) {
      conservation <- cons$ms$frequency
    } else {
      conservation <- rep(NA, length(aa_sequence))
    }
    
    if("relASA" %in% names(dssp$ms)) {
      asa <- dssp$ms$relASA
    } else {
      asa <- rep(NA, length(aa_sequence))
    }
    
    fois <- features %>% filter(type == "W or Y")
    
    for(y in regions) {
      site <- unlist(features[features[["type"]] == y,c("start", "end")], use.names = FALSE)
      
      res <- fois %>%
        mutate(overlap = start <= site[2] & end >= site[1]) %>%
        filter(overlap) %>%
        mutate(WY_cons = conservation[start],
               WY_asa = asa[start],
               WY_total = rep(nrow(.), nrow(.)),
               WY_site = paste0(aa_sequence[start], start))
      if(nrow(res) > 0) {
        res <- res %>%
          filter(WY_cons == min(WY_cons)) %>%
          dplyr::slice(1) %>%
          dplyr::select(WY_cons, WY_asa, WY_total, WY_site)
      }
      
      if(nrow(res) == 1) {
        features[features[["type"]] == y, c("WY_cons", "WY_asa", "WY_total", "WY_site")] <- res
      }
      
      if(grepl("^dbr_", y)) {
        db_inds <- as.integer(strsplit(strsplit(y, "_")[[1]][2], "-")[[1]])
        
        features[features[["type"]] == y, "Dibasic_cons"] <- mean(conservation[db_inds], na.rm = TRUE)
        features[features[["type"]] == y, "Dibasic_asa"] <-  mean(asa[db_inds], na.rm = TRUE)
      }
    }
    
  }
  
  return(features)
  
}

#secretome <- secretome %>%
#  mutate(features = pmap(list(features, cons, dssp, sequence_uni), .f = add_site_stats))


saveRDS(secretome, paste0(s_localDir, "/processed/secretome_6.rds"))


