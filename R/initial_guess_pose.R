## ---------------------------------------------------------------------------
## Initial-guess pose generation for receptor + ligand complexes.
##
## Places a ligand monomer PDB next to a receptor PDB so that either the N- or
## C-terminus sits above the receptor's extracellular (EC) face, with the bulk
## of the ligand extending outward and no inter-chain heavy-atom clashes.
##
## Designed to produce PDBs consumable by the AlphaPulldown "initial guess"
## extension (see ~/Desktop/ap_initial_guess/).  This file is standalone — it
## does not call compute_RLdists() or any of the metric-processing machinery.
## The caller supplies vectors of receptor residue numbers at the EC and IC
## faces directly (typically the EC/IC ends of TM helices).
##
## Requires: a ligand monomer PDB (fold first if you only have sequence).
## ---------------------------------------------------------------------------


# ---- geometry helpers ----------------------------------------------------

.unit <- function(v) {
  n <- sqrt(sum(v * v))
  if (n < 1e-9) stop("zero-length vector")
  v / n
}

# Rodrigues' formula: build a 3x3 rotation taking unit vector a -> unit vector b.
.rot_a_to_b <- function(a, b) {
  a <- .unit(a); b <- .unit(b)
  d <- sum(a * b)
  if (d >  1 - 1e-9) return(diag(3))                       # already aligned
  if (d < -1 + 1e-9) {                                     # antiparallel
    ortho <- if (abs(a[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
    axis  <- .unit(ortho - sum(ortho * a) * a)
    K <- matrix(c( 0, -axis[3],  axis[2],
                   axis[3],  0, -axis[1],
                  -axis[2],  axis[1],  0), 3, 3, byrow = TRUE)
    return(diag(3) + 2 * K %*% K)
  }
  axis  <- .unit(c(a[2]*b[3] - a[3]*b[2],
                   a[3]*b[1] - a[1]*b[3],
                   a[1]*b[2] - a[2]*b[1]))
  angle <- acos(d)
  K <- matrix(c( 0, -axis[3],  axis[2],
                 axis[3],  0, -axis[1],
                -axis[2],  axis[1],  0), 3, 3, byrow = TRUE)
  diag(3) + sin(angle) * K + (1 - cos(angle)) * (K %*% K)
}

# Rotate about an arbitrary unit axis through the origin.
.rot_axis <- function(axis, angle_rad) {
  axis <- .unit(axis)
  K <- matrix(c( 0, -axis[3],  axis[2],
                 axis[3],  0, -axis[1],
                -axis[2],  axis[1],  0), 3, 3, byrow = TRUE)
  diag(3) + sin(angle_rad) * K + (1 - cos(angle_rad)) * (K %*% K)
}

# Min pairwise distance + count of pairs below cutoff between two coord matrices.
.clash_stats <- function(xyz_a, xyz_b, cutoff) {
  # xyz_a: (Na, 3); xyz_b: (Nb, 3).  Vectorised, no loops.
  da2 <- rowSums(xyz_a * xyz_a)
  db2 <- rowSums(xyz_b * xyz_b)
  cross <- xyz_a %*% t(xyz_b)
  d2 <- outer(da2, db2, `+`) - 2 * cross
  d2[d2 < 0] <- 0
  list(min_dist = sqrt(min(d2)),
       n_clashes = sum(d2 < cutoff * cutoff))
}


# ---- core --------------------------------------------------------------

#' Position a ligand near the extracellular face of a receptor.
#'
#' Rigid-body places \code{ligand_pdb} relative to \code{receptor_pdb} so that
#' the chosen terminus (N or C) sits above the receptor's EC face and the bulk
#' of the ligand extends outward into the extracellular space.  Searches over
#' axial spins and clearance values to find a clash-free placement.
#'
#' @param receptor_pdb  Path to receptor PDB.
#' @param ligand_pdb    Path to ligand monomer PDB.
#' @param ec_residues   Integer vector of receptor residue numbers at the EC
#'   face (typically EC ends of TM helices).
#' @param ic_residues   Integer vector of receptor residue numbers at the IC
#'   face (typically IC ends of TM helices).
#' @param terminus      "N" or "C" — which ligand end sits near the receptor.
#' @param receptor_chain Output chain ID for the receptor (default "A").
#' @param ligand_chain   Output chain ID for the ligand   (default "B").
#' @param clearance     Initial clearance above the EC centroid, in Angstroms.
#' @param max_clearance Maximum clearance before giving up.
#' @param clash_cutoff  Heavy-atom clash distance, in Angstroms.
#' @param spin_steps    Number of axial-spin orientations to try at each
#'   clearance (evenly spaced around the membrane normal).
#' @param output_pdb    Path to write the combined PDB.  If \code{NULL}, no
#'   file is written; the returned bio3d object is the only output.
#'
#' @return A list with \code{$pdb} (bio3d combined receptor+ligand),
#'   \code{$diagnostics} (clearance/spin/min_dist/n_clashes used), and
#'   \code{$success} (TRUE if a clash-free placement was found).
#' @export
position_ligand_initial_guess <- function(
    receptor_pdb,
    ligand_pdb,
    ec_residues,
    ic_residues,
    terminus       = c("N", "C"),
    receptor_chain = "A",
    ligand_chain   = "B",
    clearance      = 5,
    max_clearance  = 25,
    clash_cutoff   = 2.0,
    spin_steps     = 12,
    output_pdb     = NULL) {

  terminus <- match.arg(terminus)

  rec <- bio3d::read.pdb(receptor_pdb, verbose = FALSE)
  lig <- bio3d::read.pdb(ligand_pdb,   verbose = FALSE)

  ## --- receptor anchor geometry ---
  ca_rec <- rec$atom[rec$atom$elety == "CA", ]
  ec_xyz <- as.matrix(ca_rec[ca_rec$resno %in% ec_residues, c("x", "y", "z")])
  ic_xyz <- as.matrix(ca_rec[ca_rec$resno %in% ic_residues, c("x", "y", "z")])
  if (nrow(ec_xyz) == 0) stop("no receptor CA atoms matched ec_residues")
  if (nrow(ic_xyz) == 0) stop("no receptor CA atoms matched ic_residues")

  ec_mid <- colMeans(ec_xyz)
  ic_mid <- colMeans(ic_xyz)
  axis   <- .unit(ic_mid - ec_mid)   # points IC; outward (EC) = -axis
  outward <- -axis

  ## --- ligand: identify chosen terminus and centre of mass ---
  ca_lig <- lig$atom[lig$atom$elety == "CA", ]
  if (nrow(ca_lig) < 2) stop("ligand has fewer than 2 CA atoms")
  ca_lig <- ca_lig[order(ca_lig$resno), ]
  term_xyz <- as.numeric(ca_lig[if (terminus == "N") 1L else nrow(ca_lig),
                                c("x", "y", "z")])
  lig_com  <- colMeans(as.matrix(ca_lig[, c("x", "y", "z")]))

  ## --- alignment: rotate ligand so (terminus -> COM) aligns with outward ---
  lig_xyz <- matrix(lig$xyz, ncol = 3, byrow = TRUE)
  cur_dir <- .unit(lig_com - term_xyz)
  R0      <- .rot_a_to_b(cur_dir, outward)
  centred <- sweep(lig_xyz, 2, term_xyz, "-")   # terminus at origin
  base_rotated <- centred %*% t(R0)

  rec_heavy <- as.matrix(rec$atom[rec$atom$elety != "H" &
                                    !startsWith(rec$atom$elety, "H"),
                                  c("x", "y", "z")])

  ## --- search: increasing clearance, axial spin around outward normal ---
  best <- list(min_dist = -Inf, n_clashes = NA_integer_,
               clearance = NA_real_, spin = NA_real_, xyz = NULL)

  step_size <- 1
  cls <- clearance
  while (cls <= max_clearance) {
    anchor <- ec_mid + cls * outward
    for (k in 0:(spin_steps - 1)) {
      angle  <- 2 * pi * k / spin_steps
      Rspin  <- .rot_axis(outward, angle)
      placed <- base_rotated %*% t(Rspin)
      placed <- sweep(placed, 2, anchor, "+")
      stats  <- .clash_stats(rec_heavy, placed, clash_cutoff)
      if (stats$n_clashes == 0) {
        best <- list(min_dist = stats$min_dist, n_clashes = 0L,
                     clearance = cls, spin = angle, xyz = placed)
        break
      }
      if (stats$min_dist > best$min_dist) {
        best <- list(min_dist = stats$min_dist, n_clashes = stats$n_clashes,
                     clearance = cls, spin = angle, xyz = placed)
      }
    }
    if (!is.na(best$n_clashes) && best$n_clashes == 0) break
    cls <- cls + step_size
  }

  ## --- assemble combined PDB ---
  # best$xyz rows are in atom-record order, same order as lig$atom.
  lig_out <- lig
  lig_out$atom$x <- best$xyz[, 1]
  lig_out$atom$y <- best$xyz[, 2]
  lig_out$atom$z <- best$xyz[, 3]
  lig_out$xyz <- matrix(t(best$xyz), nrow = 1)

  rec_atoms <- rec$atom
  lig_atoms <- lig_out$atom
  rec_atoms$chain <- receptor_chain
  lig_atoms$chain <- ligand_chain
  lig_atoms$eleno <- lig_atoms$eleno + max(rec_atoms$eleno)

  combined_atoms <- rbind(rec_atoms, lig_atoms)
  combined_xyz   <- as.numeric(t(as.matrix(combined_atoms[, c("x", "y", "z")])))

  combined_pdb <- list(atom = combined_atoms,
                       xyz  = matrix(combined_xyz, nrow = 1),
                       helix = NULL, sheet = NULL, seqres = NULL,
                       calpha = combined_atoms$elety == "CA",
                       call = match.call())
  class(combined_pdb) <- c("pdb", "sse")

  if (!is.null(output_pdb)) {
    bio3d::write.pdb(pdb = combined_pdb, file = output_pdb)
  }

  list(
    pdb         = combined_pdb,
    success     = isTRUE(best$n_clashes == 0),
    diagnostics = tibble::tibble(
      terminus    = terminus,
      clearance   = best$clearance,
      spin_rad    = best$spin,
      min_dist    = best$min_dist,
      n_clashes   = best$n_clashes
    )
  )
}


# ---- batch -------------------------------------------------------------

#' Batch-position a set of receptor/ligand pairs.
#'
#' @param jobs A data frame with columns:
#'   \itemize{
#'     \item \code{pair_id}        — used for output filename (\code{<pair_id>.pdb}).
#'     \item \code{receptor_pdb}   — path to receptor PDB.
#'     \item \code{ligand_pdb}     — path to ligand monomer PDB.
#'     \item \code{ec_residues}    — list-column of integer vectors,
#'                                    OR character (comma-separated).
#'     \item \code{ic_residues}    — same format as \code{ec_residues}.
#'     \item \code{terminus}       — "N" or "C".
#'   }
#' @param output_dir Directory to write PDBs into (created if missing).
#' @param ...        Forwarded to \code{position_ligand_initial_guess()}
#'                   (clearance, clash_cutoff, spin_steps, etc.).
#'
#' @return A tibble joining \code{jobs} with per-row diagnostics and an
#'   \code{output_pdb} column.
#' @export
position_ligands_initial_guess_batch <- function(jobs, output_dir, ...) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  parse_resids <- function(x) {
    if (is.list(x))      return(as.integer(x[[1]]))
    if (is.numeric(x))   return(as.integer(x))
    as.integer(strsplit(as.character(x), "[, ]+")[[1]])
  }

  rows <- lapply(seq_len(nrow(jobs)), function(i) {
    job <- jobs[i, ]
    out <- file.path(output_dir, paste0(job$pair_id, ".pdb"))
    res <- tryCatch(
      position_ligand_initial_guess(
        receptor_pdb = job$receptor_pdb,
        ligand_pdb   = job$ligand_pdb,
        ec_residues  = parse_resids(job$ec_residues),
        ic_residues  = parse_resids(job$ic_residues),
        terminus     = as.character(job$terminus),
        output_pdb   = out,
        ...
      ),
      error = function(e) list(success = FALSE,
                               diagnostics = tibble::tibble(error = conditionMessage(e)))
    )
    tibble::tibble(
      pair_id     = job$pair_id,
      output_pdb  = if (isTRUE(res$success)) out else NA_character_,
      success     = isTRUE(res$success),
      !!!res$diagnostics
    )
  })

  do.call(dplyr::bind_rows, rows)
}
