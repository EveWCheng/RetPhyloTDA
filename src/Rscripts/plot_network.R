.libPaths("~/R/library")
library(ape)
library(SiPhyNetwork)
library(yaml)

# ── Config ────────────────────────────────────────────────────────────────────
# First CLI arg is the path to a config yaml (input_dir / output_dir). Any
# further args are sim indices to restrict the run to; if none are given,
# every simN_nodes.csv/simN_edges.csv pair found in input_dir is plotted.
# Usage: Rscript plot_network.R <config.yaml> [sim_indices...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_network.R <config.yaml> [sim_indices...]")
}
config   <- yaml::read_yaml(args[1])
IN_DIR   <- config$input_dir
PLOT_DIR <- config$output_dir
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

sim_args <- args[-1]
if (length(sim_args) > 0) {
  sim_indices <- as.integer(sim_args)
} else {
  node_files  <- list.files(IN_DIR, pattern = "^sim[0-9]+_nodes\\.csv$")
  sim_indices <- sort(as.integer(gsub("sim([0-9]+)_nodes\\.csv", "\\1", node_files)))
}

# ── Build an ape/evonet phylo object from one sim's nodes/edges CSVs ──────────

build_phy <- function(nodes, edges) {
  tree_edges <- edges[edges$edge_type == "tree", ]
  ret_edges  <- edges[edges$edge_type == "reticulation", ]

  # Tips: nodes with no outgoing tree edges (leaf, extinct, and hyb_source
  # nodes all qualify - none of them have tree-edge children).
  tip_ids <- nodes$id[!(nodes$id %in% tree_edges$from)]
  # Root: node with no incoming edges (not a child in any edge)
  root_id <- nodes$id[!(nodes$id %in% c(tree_edges$to, ret_edges$to))]
  # Internal: everything else (put root first so it maps to ntips + 1)
  internal_ids <- c(root_id, setdiff(nodes$id[nodes$id %in% tree_edges$from], root_id))

  ntips <- length(tip_ids)
  Nnode <- length(internal_ids)

  # Mapping: sim node ID -> ape node ID
  id_map <- c(setNames(seq_len(ntips),              as.character(tip_ids)),
              setNames(seq(ntips + 1, ntips + Nnode), as.character(internal_ids)))

  edge_mat <- matrix(c(id_map[as.character(tree_edges$from)],
                       id_map[as.character(tree_edges$to)]),
                     ncol = 2)

  tip_labels <- nodes$label[match(tip_ids, nodes$id)]

  phy <- structure(
    list(edge        = edge_mat,
         edge.length = tree_edges$length,
         tip.label   = tip_labels,
         Nnode       = Nnode),
    class = "phylo"
  )
  phy <- reorder(phy)

  if (nrow(ret_edges) > 0) {
    phy$reticulation <- matrix(
      c(id_map[as.character(ret_edges$from)],
        id_map[as.character(ret_edges$to)]),
      ncol = 2)
    phy$reticulation.length <- ret_edges$length
    class(phy) <- c("evonet", "phylo")
  }

  attr(phy, "internal_ids") <- internal_ids
  phy
}

# ── Read one sim's CSVs, plot it, and write the PDF ───────────────────────────

plot_sim <- function(sim_index) {
  nodes_file <- file.path(IN_DIR, paste0("sim", sim_index, "_nodes.csv"))
  edges_file <- file.path(IN_DIR, paste0("sim", sim_index, "_edges.csv"))

  nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
  edges <- read.csv(edges_file, stringsAsFactors = FALSE)

  phy <- build_phy(nodes, edges)
  internal_ids <- attr(phy, "internal_ids")

  out_file <- file.path(PLOT_DIR, paste0("sim", sim_index, "_network.pdf"))
  pdf(out_file)
  on.exit(dev.off(), add = TRUE)

  if (length(phy$tip.label) < 2) {
    # ape::plot.phylo/nodelabels require >= 2 tips to set up a plot region
    plot.new()
    title(main = paste("BDH Network - sim", sim_index))
    text(0.5, 0.5, paste("Only", length(phy$tip.label), "tip - nothing to plot"))
  } else {
    plot(phy, main = paste("BDH Network - sim", sim_index))
    internal_labels <- nodes$label[match(internal_ids, nodes$id)]
    ape::nodelabels(text = internal_labels, cex = 0.6, frame = "none", col = "blue")
  }

  cat("Plot written to", out_file, "\n")
}

# ── Run ────────────────────────────────────────────────────────────────────────

for (sim_index in sim_indices) {
  result <- tryCatch({
    plot_sim(sim_index)
    TRUE
  }, error = function(e) {
    cat("Skipping sim", sim_index, "- error:", conditionMessage(e), "\n")
    FALSE
  })
}
