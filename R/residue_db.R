





set_db_path <- function(path) {
  config_file <- file.path(tools::R_user_dir("ligandFinder", which = "config"), "config.txt")
  dir.create(dirname(config_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(path, config_file)
}



ligandFinder::get_db_path <- function() {
  config_file <- file.path(tools::R_user_dir("ligandFinder", which = "config"), "config.txt")

  if (file.exists(config_file)) {
    return(readLines(config_file, warn = FALSE)[1])
  }

  return(tools::R_user_dir("ligandFinder", which = "data"))  # Default path
}


download_residue_db <- function(dest_dir = ligandFinder:::ligandFinder::get_db_path()) {

  db_url <- "https://stacks.stanford.edu/file/sc075gg6264/residue_db.tar.gz"

  db_path <- fs::path_expand(paste0(dest_dir, "/residue_db.tar.gz"))

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }

  if (!file.exists(db_path)) {
    message("Downloading database...")
    download.file(db_url, db_path, mode = "wb")
    untar(db_path, exdir = fs::path_expand(dest_dir))
    fs::file_delete(db_path)
    message("Database downloaded to: ", paste0(dest_dir, "/residue_db"))
  } else {
    message("Database already exists at: ", db_path)
  }

  return(db_path)
}


