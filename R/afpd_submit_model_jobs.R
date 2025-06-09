


download_roi_data <- function(url,
                              dest_dir = get_db_path()) {

  downloaded_path <- paste0(dest_dir, "/", stringr::str_extract(url, "[^/]+$"))


  if (!file.exists(downloaded_path)) {
    message("Downloading data...")
    download.file(url, downloaded_path, mode = "wb")
    message("data downloaded to: ", downloaded_path)
  } else {
    message("data already exists at: ", downloaded_path)
  }

  return(invisible(downloaded_path))
}



submit_model_jobs <- function(script = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh",
                              input_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/mxrun",
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
          system(paste("sbatch",
                       paste0("--job-name=", sub(".txt$", "", jobs[ind])),
                       script,
                       paste0(input_path, "/", jobs[ind]),
                       output_path))
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
