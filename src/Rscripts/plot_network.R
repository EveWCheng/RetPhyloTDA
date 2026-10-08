.libPaths("~/R/library")
library(ape)
library(SiPhyNetwork)
library(yaml)

# ── Config ────────────────────────────────────────────────────────────────────
# First CLI arg is the path to a config yaml (input_dir / output_dir). Any
# further args are sim indices to restrict the run to; if none are given,
# every sim with a simN_nodes.csv or simN_filtered_nodes.csv in input_dir is
# plotted. Relative paths in the config are resolved against the config file's
# own directory, so the script can be run from any working directory.
# Usage: Rscript plot_network.R <config.yaml> [sim_indices...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_network.R <config.yaml> [sim_indices...]")
}
config     <- yaml::read_yaml(args[1])
CONFIG_DIR <- dirname(normalizePath(args[1], mustWork = TRUE))
resolve_path <- function(p) {
  if (grepl("^(/|~)", p)) path.expand(p) else file.path(CONFIG_DIR, p)
}
IN_DIR   <- normalizePath(resolve_path(config$input_dir), mustWork = TRUE)
PLOT_DIR <- resolve_path(config$output_dir)
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)
PLOT_DIR <- normalizePath(PLOT_DIR)
# Optional: main_result_analysis.py's simN_analysis.txt reports, used to colour detected reticulations
REPORT_DIR <- if (!is.null(config$report_dir)) resolve_path(config$report_dir) else NULL

sim_args <- args[-1]
if (length(sim_args) > 0) {
  sim_indices <- as.integer(sim_args)
} else {
  # sim_main.py only exports the filtered CSVs; tree_main.py exports both.
  node_files  <- list.files(IN_DIR, pattern = "^sim[0-9]+_(filtered_)?nodes\\.csv$")
  sim_indices <- sort(unique(as.integer(gsub("^sim([0-9]+)_.*$", "\\1", node_files))))
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

# ── Reticulations detected per main_result_analysis.py's report ──────────────

# Returns "from--to" keys of reticulations listed as on / not on a cycle-edge path, or NULL
# when there is no report. Reads whatever top-k the report was written with.
read_detection <- function(sim_index) {
  if (is.null(REPORT_DIR)) return(NULL)
  report_file <- file.path(REPORT_DIR, paste0("sim", sim_index, "_analysis.txt"))
  if (!file.exists(report_file)) return(NULL)
  lines <- readLines(report_file)
  on_start  <- grep("^reticulate edges on at least one top-[0-9]+ cycle-edge path", lines)
  off_start <- grep("^reticulate edges on no top-[0-9]+ cycle-edge path", lines)
  if (length(on_start) != 1 || length(off_start) != 1) return(NULL)
  edge_keys <- function(block) {
    hits <- regmatches(block, regexec("^  (\\S+) -- (\\S+)  inher_weight=", block))
    unlist(lapply(hits, function(m) if (length(m) == 3) paste(m[2], m[3], sep = "--")))
  }
  off_end <- off_start + which(!grepl("^  ", lines[-seq_len(off_start)]))[1] - 1
  if (is.na(off_end)) off_end <- length(lines)
  list(detected = edge_keys(lines[(on_start + 1):(off_start - 1)]),
       missed   = edge_keys(lines[seq(off_start + 1, length.out = off_end - off_start)]))
}

# ── Read one sim's CSVs, plot it, and write the PDF ───────────────────────────

# kind = "" plots simN_nodes.csv/simN_edges.csv (raw phy.G, export_csv);
# kind = "filtered_" plots simN_filtered_nodes.csv/simN_filtered_edges.csv
# (the "no_hyb_nodes"-collapsed graph gene-tree enumeration actually runs on,
# export_filtered) -- same node IDs as the raw graph, just fewer nodes, so the
# two PDFs stay directly comparable.
plot_sim <- function(sim_index, kind = "") {
  nodes_file <- file.path(IN_DIR, paste0("sim", sim_index, "_", kind, "nodes.csv"))
  edges_file <- file.path(IN_DIR, paste0("sim", sim_index, "_", kind, "edges.csv"))

  nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
  edges <- read.csv(edges_file, stringsAsFactors = FALSE)

  phy <- build_phy(nodes, edges)
  internal_ids <- attr(phy, "internal_ids")

  out_file <- file.path(PLOT_DIR, paste0("sim", sim_index, "_", kind, "network.pdf"))
  pdf(out_file)
  on.exit(dev.off(), add = TRUE)

  title_label <- paste0("BDH Network - ", if (kind == "") "raw" else "filtered", " - sim", sim_index)

  if (length(phy$tip.label) < 2) {
    # ape::plot.phylo/nodelabels require >= 2 tips to set up a plot region
    plot.new()
    title(main = title_label)
    text(0.5, 0.5, paste("Only", length(phy$tip.label), "tip - nothing to plot"))
  } else {
    # reticulations: red = on a cycle-edge path in the analysis report, grey dashed = missed,
    # blue = no report (same order as build_phy's phy$reticulation rows)
    ret_edges <- edges[edges$edge_type == "reticulation", ]
    detection <- if (kind == "filtered_") read_detection(sim_index) else NULL
    ret_col <- "blue"; ret_lty <- 1; ret_lwd <- 1
    if (!is.null(detection) && nrow(ret_edges) > 0) {
      fwd <- paste(ret_edges$from, ret_edges$to, sep = "--")
      rev <- paste(ret_edges$to, ret_edges$from, sep = "--")
      hit  <- fwd %in% detection$detected | rev %in% detection$detected
      miss <- fwd %in% detection$missed   | rev %in% detection$missed
      ret_col <- ifelse(hit, "red", ifelse(miss, "grey50", "blue"))
      ret_lty <- ifelse(hit, 1, 2)
      ret_lwd <- 1
      title_label <- paste0(title_label, " - ", sum(hit), "/", nrow(ret_edges), " retic detected")
    }
    plot(phy, main = title_label, col = ret_col, lty = ret_lty, lwd = ret_lwd, alpha = 0.9)
    internal_labels <- nodes$label[match(internal_ids, nodes$id)]
    ape::nodelabels(text = internal_labels, cex = 0.6, frame = "none", col = "blue")
    # edge lengths (2 dp); zero-length edges are left unlabelled to keep the plot readable
    tree_len <- phy$edge.length
    ape::edgelabels(text = ifelse(tree_len > 0, formatC(tree_len, format = "f", digits = 2), ""),
                    cex = 0.5, frame = "none", adj = c(0.5, -0.4), col = "darkgreen")
    if (!is.null(phy$reticulation)) {
      # ape has no labeller for reticulation lines, so place each length at its line's midpoint
      lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
      mid_x <- (lp$xx[phy$reticulation[, 1]] + lp$xx[phy$reticulation[, 2]]) / 2
      mid_y <- (lp$yy[phy$reticulation[, 1]] + lp$yy[phy$reticulation[, 2]]) / 2
      text(mid_x, mid_y, formatC(phy$reticulation.length, format = "f", digits = 2),
           cex = 0.5, col = ifelse(is.character(ret_col) & ret_col == "red", "red", "grey30"))
    }
    if (!is.null(detection) && nrow(ret_edges) > 0) {
      legend("topleft", legend = c("retic detected", "retic missed"), col = c("red", "grey50"),
             lty = c(1, 2), lwd = c(1, 1), bty = "n", cex = 0.7)
    }
  }

  cat("Plot written to", out_file, "\n")
}

# ── Run ────────────────────────────────────────────────────────────────────────

for (sim_index in sim_indices) {
  for (kind in c("", "filtered_")) {
    # Not every pipeline exports both kinds, so skip silently when absent
    if (!file.exists(file.path(IN_DIR, paste0("sim", sim_index, "_", kind, "nodes.csv")))) next
    result <- tryCatch({
      plot_sim(sim_index, kind)
      TRUE
    }, error = function(e) {
      cat("Skipping sim", sim_index, "kind='", kind, "' - error:", conditionMessage(e), "\n")
      FALSE
    })
  }
}
