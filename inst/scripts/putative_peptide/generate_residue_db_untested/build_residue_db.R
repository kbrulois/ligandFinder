

library(tidyverse)
library(ligandFinder)
residue_db_scripts <- system.file("/scripts/putative_peptide/generate_residue_db_untested", package = "ligandFinder")

sorce <- function(x) {source(fs::path(residue_db_scripts, x))}

sorce("2_parse_aminode.R")

sorce("2.1_parse_aminode_new.R")

sorce("3_integrate_aminode.R")

sorce("4_integrate_dssp.R")

sorce("5_integrate_afm.R")

sorce("6_integrate_peptides.R")

sorce("7_define_phs.R")

sorce("8_define_ec.R")


