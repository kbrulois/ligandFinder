

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)


submit_model_jobs <- function(script = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh",
                        protein_input_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/mxrun",
                        output_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/mxrun",
                        jobs = unique(to_run$group),
                        start_from = 1,
                        max_jobs = 16,
                        sleep_time = 60) {

  jobs <- jobs[start_from:length(jobs)]
  i <- 1
  total_jobs <- length(jobs)

  repeat{
        current_jobs <- system("squeue -u $USER | grep -c 'gpu'", intern = TRUE)
        if(stringr::str_detect(current_jobs, "^[0-9]+$")) {
          current_jobs <- as.numeric(current_jobs)
        } else {
          current_jobs <- max_jobs
        }
          if (current_jobs < max_jobs & current_jobs >= 0) {
              message(current_jobs, " current jobs pending or running")
              num_jobs_to_add <- max_jobs - current_jobs
              for(x in 1:num_jobs_to_add) {
                  ind <- i + x - 1
                  if(ind <= total_jobs) {
                  system(paste("sbatch", script, paste0(protein_input_path, "/", jobs[ind]), output_path))
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

            start_time <- Sys.time()
            #do fake calculation so your job doesn't get killed
            for (num in 1:1e8) {
              sqrt(num)
              if (difftime(Sys.time(), start_time, units = "secs") > sleep_time) break
            }


          }
      }
}

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

gpcr_sub <- gpcr_list %>%
              filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
              filter(`ecb: Exclude due to N term >160AA` != "Exclude due to N term >160AA") %>%
              filter(map_lgl(`bw: full_table`, ~nrow(.) > 0))

ligand_list <- tibble(uniprot_name = c("GP15L", "GP15L", "CXL17", "CCL25", "RARR2", "RARR2", "RARR2", "RARR2"),
                      start = c(25, 71, 64, 24, 21, 21, 21,21),
                      end = c(81, 81, 119, 150, 163, 157, 158, 156))

ligand_list <- tibble(end = 88:107,
                      uniprot_name = "CXL14",
                      start = 35) %>%

ligand_list <- tibble(uniprot_name = c("CART", "CART", "CART", "CART", "CQ067", "CQ067", "SDF1", "SDF1", "CO061", "CO061", "CO061", "RARR2", "CCL27", "DMKN", "DMKN", "DMKN"),
                      start = c(28, 64, 76, 93, 19, 38, 71, 71, 29, 31, 45, 29, 77, 460, 355, 395),
                      end = c(61, 73, 90, 116, 35, 90, 82, 88, 42, 42, 157, 45, 88, 476, 376, 408),
                      run = "mxrun")

ligand_list <- tibble(uniprot_name = c("LUZP2", "LUZP2", "LUZP2", "LUZP2"),
                      start = c(157, 123, 22, 327),
                      end = c(167, 154, 62, 346),
                      run = "luzp2")



ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))

ligand_list <- ligand_list %>%
                  mutate(model = paste0(uniprot_name, ",", start, "-", end)) %>%
                  filter(ecb_cull == "y" & database %in% c("both", "gpcrdb"))

ligand_list <- left_join(ligand_list, gpcr_list %>% dplyr::rename(receptor = uniprot_name), by = "receptor") %>%
                  mutate(known_model = paste0(model_name, ";", model))


to_run <- expand.grid(ligand = ligand_list[["model"]],
                      receptor = gpcr_sub[["model_name"]],
                      stringsAsFactors = FALSE) %>%
          as_tibble %>%
          mutate(model = paste(receptor, ligand, sep = ";"))


to_run <- to_run %>%
              mutate(known = if_else(model %in% ligand_list$known_model, "known", "unknown")) %>%
              arrange(known)

group_size <- 48

to_run <- to_run %>%
    mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))


to_run %>%
      group_by(group) %>%
      group_walk(~ write.table(.x$model, file = .y$group,
                               row.names = FALSE, col.names = FALSE, quote = FALSE))

job_script <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh"

submit_model_jobs(script = job_script,
                  protein_input_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/benchmarking",
                  output_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking",
                  jobs = unique(to_run$group),
                  start_from = 23,
                  max_jobs = 16)






comp_jobs <- tibble(afpd_file_name = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking"))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

comp_jobs <- comp_jobs %>%
  mutate(parse_proteins(afpd_file_name)) %>%
  mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
  mutate(model = paste0(p1_id,
                        if_else(p1_range == "", "", paste0(",", p1_range)),
                        ";",
                        p2_id,
                        if_else(p2_range == "", "", paste0(",", p2_range))))

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


submit_jobs(script = job_script,
            jobs = unique(to_still_run$group),
            start_from = 1,
            max_jobs = 16)













