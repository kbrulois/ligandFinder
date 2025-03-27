

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)

gpcr_list <- readRDS("/oak/stanford/groups/ebutcher/deorphan-AI-ze/uniprot/gpcrdb_receptor_list_KB_250304_final.rds")

gpcr_sub <- gpcr_list %>%
  filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  filter(`ecb: Exclude due to N term >160AA` != "Exclude due to N term >160AA") %>%
  #mutate(gpcr_priority = if_else(gene_name_primary %in% gpcr_ord[["V1"]], "group1", "group2")) %>%
  #group_by(gpcr_priority) %>%
  #mutate(order = match(gene_name_primary, gpcr_ord[["V1"]])) %>%
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  #arrange(gpcr_priority, order) %>%
  select(uniprot_name, 
         `bw: indices 1.50`, 
         `bw: indices 7.50`, 
         contains("ecb: Length"),
         signal_peptide,
         last_AA) %>%
  mutate(model = if_else(is.na(signal_peptide), 
                         uniprot_name, 
                         paste0(uniprot_name, ",", signal_peptide, "-", last_AA)))


ligand_list <- tibble(uniprot_name = c("GP15L", "GP15L", "CXL17", "CCL25", "RARR2", "RARR2", "RARR2", "RARR2"),
                      start = c(25, 71, 64, 24, 21, 21, 21,21),
                      end = c(81, 81, 119, 150, 163, 157, 158, 156)) %>%
  mutate(model = paste0(uniprot_name, ",", start, "-", end))


group_size <- 32

to_run <- expand.grid(ligand = ligand_list[["model"]], 
                      receptor = gpcr_sub[["model"]], 
                      stringsAsFactors = FALSE) %>%
  as_tibble %>%
  mutate(model = paste(receptor, ligand, sep = ";")) %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x$model, file = .y$group, 
                           row.names = FALSE, col.names = FALSE, quote = FALSE))



submit_jobs <- function(script = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown_BRINP2.sh", 
                        jobs = unique(to_run$group), 
                        max_jobs = 20, 
                        sleep_time = 60) {
  
  i <- 1
  total_jobs <- length(jobs)
  
  on.exit(return(i))
  
  repeat{
    current_jobs <- as.numeric(system("squeue -u $USER | wc -l", intern = TRUE)) - 1
    
    if (current_jobs < max_jobs) {
      
      message(current_jobs, " current jobs pending or running")
      
      num_jobs_to_add <- max_jobs - current_jobs
      for(x in 1:num_jobs_to_add) {
        
        ind <- i + x - 1
        if(ind < total_jobs) {
          system(paste("sbatch", script, jobs[ind]))
          message("Submitted a new job: ", jobs[ind])
        } else {
          break
        }
      }
      i <- i + num_jobs_to_add
      message("i: ", i)
      if(i > total_jobs) {
        break
      }
      message("submitted ", (i - 1), " of ", total_jobs, " total jobs")
      message("Jobs are currently maxed. Waiting...")
      Sys.sleep(30)  
    } else {
      Sys.sleep(sleep_time)  
    }
  }
  
}

submit_jobs()






comp_jobs <- tibble(afpd_file_name = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/pocForGrant"))

#comp_jobs <- tibble(files = list.files("~/oak/deorphan-AI-ze/models/pocForGrant"))

id_map <- readRDS("/oak/stanford/groups/ebutcher/deorphan-AI-ze/uniprot/id_mapping.rds")

comp_jobs <- comp_jobs %>%
  mutate(parse_proteins(files)) %>%
  mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
  mutate(model = paste0(receptor_id, 
                        if_else(receptor_range == "", "", paste0(",", receptor_range)), 
                        ";", 
                        ligand_id, 
                        if_else(ligand_range == "", "", paste0(",", ligand_range))))

sum(to_run[["model"]] %in% comp_jobs[["model"]])


to_still_run <- to_run %>%
  filter(!model %in% comp_jobs[["model"]]) %>%
  mutate(parse_proteins(model, delim_proteins = ";", delim_ranges = ","))

left_join(to_still_run, id_map %>% rename(receptor_id = `Entry Name`), by = "receptor_id") %>%
  select(receptor_id, Length) %>%
  print(n = 1000)




to_still_run <- to_still_run %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n())) 

to_still_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x$model, file = .y$group, 
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


submit_jobs <- function(script = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh", 
                        jobs = unique(to_still_run$group), 
                        max_jobs = 20, 
                        sleep_time = 60) {
  
  i <- 1
  total_jobs <- length(jobs)
  
  on.exit(return(i))
  
  repeat{
    current_jobs <- as.numeric(system("squeue -u $USER | wc -l", intern = TRUE)) - 1
    
    if (current_jobs < max_jobs) {
      
      message(current_jobs, " current jobs pending or running")
      
      num_jobs_to_add <- max_jobs - current_jobs
      for(x in 1:num_jobs_to_add) {
        
        ind <- i + x - 1
        if(ind <= total_jobs) {
          system(paste("sbatch", script, jobs[ind]))
          message("Submitted a new job: ", jobs[ind])
        } else {
          break
        }
      }
      i <- i + num_jobs_to_add
      message("i: ", i)
      if(i > total_jobs) {
        break
      }
      message("submitted ", (i - 1), " of ", total_jobs, " total jobs")
      message("Jobs are currently maxed. Waiting...")
      Sys.sleep(30)  
    } else {
      Sys.sleep(sleep_time)  
    }
  }
  
}

submit_jobs()













