# ── Dependencies ──────────────────────────────────────────────────────────────
.libPaths("~/R/library")
library(SiPhyNetwork)

# ── Load simulation functions ──────────────────────────────────────────────────
source("sim.bdh.age.help.2_update.R")

# ── Functions ──────────────────────────────────────────────────────────────────

run_simulation <- function(age, lambda, mu, nu, hybprops,
                           hyb.inher.fxn, hyb.rate.fxn = NULL,
                           frac = 1, mrca = FALSE, complete = TRUE,
                           stochsampling = FALSE, trait.model = NULL,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  with_hyb <- sim.bdh.age.help2(
    dummy         = 1,
    age           = age,
    lambda        = lambda,
    mu            = mu,
    nu            = nu,
    hybprops      = hybprops,
    hyb.rate.fxn  = hyb.rate.fxn,
    hyb.inher.fxn = hyb.inher.fxn,
    frac          = frac,
    mrca          = mrca,
    complete      = complete,
    stochsampling = stochsampling,
    trait.model   = trait.model
  )

  with_hyb
}

print_summary <- function(result, age) {
  phy <- result$phy

  if (is.numeric(phy) && phy == 0) {
    cat("Tree went extinct before age", age, "\n")
    return(invisible(NULL))
  }

  cat("── Simulation results ──────────────────────────────\n")
  cat("Number of extant tips     :", phy$nleaves, "\n")
  cat("Number of extinct tips    :", length(phy$extinct), "\n")
  cat("Number of hybridisations  :", nrow(phy$reticulation), "\n")
}

get_node_labels <- function(result) {
  phy        <- result$phy
  node_names <- names(result$distance)
  name_table <- result$name_table
  name_table <- name_table[!is.na(name_table[, 1]), , drop = FALSE]

  ntips      <- phy$nleaves + length(phy$extinct)
  sp_counter <- 1
  labels       <- character(length(node_names))
  ape_to_label <- c()

  for (i in seq_along(node_names)) {
    ape_num <- as.numeric(name_table[name_table[, 1] == as.numeric(node_names[i]), 2])
    if (ape_num <= ntips) {
      lbl <- paste0("sp", sp_counter)
      sp_counter <- sp_counter + 1
    } else {
      lbl <- node_names[i]
    }
    labels[i] <- lbl
    ape_to_label[as.character(ape_num)] <- lbl
  }
  list(labels = labels, ape_to_label = ape_to_label)
}

write_distance_matrix <- function(result, file = "distance_matrix.txt", all_nodes = FALSE) {
  nl         <- get_node_labels(result)
  node_names <- names(result$distance)
  labels     <- nl$labels

  if (all_nodes) {
    sel_names  <- node_names
    sel_labels <- labels
  } else {
    tip_idx    <- grep("^sp", labels)
    sel_names  <- node_names[tip_idx]
    sel_labels <- labels[tip_idx]
  }

  dist_df <- do.call(rbind, lapply(result$distance[sel_names], function(row) {
    unlist(row[sel_names])
  }))
  rownames(dist_df) <- sel_labels
  colnames(dist_df) <- sel_labels

  write.table(round(dist_df, 4),
              file      = file,
              quote     = FALSE,
              sep       = "\t",
              col.names = NA)
  cat("Distance matrix written to", file, "\n")
}


descendant_tips <- function(edge_matrix, node, ntips) {
  tips  <- c()
  queue <- node
  while (length(queue) > 0) {
    cur      <- queue[1]
    queue    <- queue[-1]
    children <- edge_matrix[edge_matrix[, 1] == cur, 2]
    tips     <- c(tips,  children[children <= ntips])
    queue    <- c(queue, children[children >  ntips])
  }
  tips
}

print_hybridisation_edges <- function(result) {
  phy   <- result$phy
  ret   <- phy$reticulation

  if (is.null(ret) || nrow(ret) == 0) {
    cat("No hybridisation edges.\n")
    return(invisible(NULL))
  }

  ape_to_label <- get_node_labels(result)$ape_to_label
  ntips        <- phy$nleaves + length(phy$extinct)

  ret_len      <- phy$reticulation.length
  cat("── Hybridisation edges ─────\n")
  for (i in seq_len(nrow(ret))) {
    donor_lbl  <- ape_to_label[as.character(ret[i, 1])]
    recip_lbl  <- ape_to_label[as.character(ret[i, 2])]
    donor_desc <- descendant_tips(phy$edge, ret[i, 1], ntips)
    recip_desc <- descendant_tips(phy$edge, ret[i, 2], ntips)
    donor_sp   <- ape_to_label[as.character(donor_desc)]
    recip_sp   <- ape_to_label[as.character(recip_desc)]
    donor_sp   <- donor_sp[!is.na(donor_sp)]
    recip_sp   <- recip_sp[!is.na(recip_sp)]
    len_str    <- if (!is.null(ret_len)) sprintf("  len: %.3f", ret_len[i]) else ""
    cat(sprintf("  %s {%s}  ->  %s {%s}%s\n",
                donor_lbl, paste(donor_sp, collapse = ", "),
                recip_lbl, paste(recip_sp, collapse = ", "),
                len_str))
  }
}

plot_network <- function(result, file = "Rplot.pdf") {
  phy          <- result$phy
  ape_to_label <- get_node_labels(result)$ape_to_label

  ntips        <- phy$nleaves + length(phy$extinct)
  n_nodes      <- ntips + phy$Nnode
  internal_ids <- seq(ntips + 1, n_nodes)

  tip_labels <- phy$tip.label
  for (i in seq_len(ntips)) {
    lbl <- ape_to_label[as.character(i)]
    if (!is.na(lbl)) tip_labels[i] <- lbl
  }
  phy$tip.label <- tip_labels
  internal_labels <- ape_to_label[as.character(internal_ids)]

  pdf(file)
  plot(phy, main = "BDH Simulation")
  ape::nodelabels(text = internal_labels, cex = 0.6, frame = "none", col = "blue")
  if (!is.null(phy$edge.length)) {
    ape::edgelabels(text = round(phy$edge.length, 3), cex = 0.5, frame = "none", col = "darkred")  # distance-based
  }
  dev.off()
  cat("Plot written to", file, "\n")
}

# ── Parameters ────────────────────────────────────────────────────────────────
AGE      <- 4
LAMBDA   <- 0.5
MU       <- 0.1
NU       <- 0.05
FRAC     <- 1.0
HYBPROPS <- c(0, 0, 1)   # c(lineage generating, degenerative, neutral)

hyb.inher.fxn <- function() runif(1, 0, 1)
hyb.rate.fxn  <- NULL

# ── Run ───────────────────────────────────────────────────────────────────────
sim <- run_simulation(
  age           = AGE,
  lambda        = LAMBDA,
  mu            = MU,
  nu            = NU,
  hybprops      = HYBPROPS,
  hyb.inher.fxn = hyb.inher.fxn,
  hyb.rate.fxn  = hyb.rate.fxn,
  frac          = FRAC,
  seed          = 50
)

OUT_DIR <- "../../Outputs/phylo_outputs"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

print_summary(sim, AGE)
print_hybridisation_edges(sim)
write_distance_matrix(sim, file.path(OUT_DIR, "distance_matrix.txt"))
write_distance_matrix(sim, file.path(OUT_DIR, "distance_matrix_all_nodes.txt"), all_nodes = TRUE)
plot_network(sim, file.path(OUT_DIR, "Rplot.pdf"))
