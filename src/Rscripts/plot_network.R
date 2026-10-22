.libPaths("~/R/library")
library(ape)
library(SiPhyNetwork)
library(yaml)

# ── Config ────────────────────────────────────────────────────────────────────
# First CLI arg is the path to a config yaml (input_dir / output_dir). Any
# further args are sim indices to restrict the run to; if none are given,
# every simN/ folder in input_dir (holding simN_nodes.csv and/or simN_filtered_nodes.csv) is
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
# Optional: main_result_analysis.py's detection_report.txt, used to colour detected reticulations
REPORT_DIR <- if (!is.null(config$report_dir)) resolve_path(config$report_dir) else NULL

# path to one of a sim's CSVs: input_dir/simN/simN_<kind><what>.csv
sim_csv_file <- function(sim_index, kind, what) {
  file.path(IN_DIR, paste0("sim", sim_index), paste0("sim", sim_index, "_", kind, what, ".csv"))
}

sim_args <- args[-1]
if (length(sim_args) > 0) {
  sim_indices <- as.integer(sim_args)
} else {
  # one simN/ folder per sim; sim_main.py only exports the filtered CSVs, tree_main.py exports both.
  sim_dirs    <- list.files(IN_DIR, pattern = "^sim[0-9]+$")
  sim_indices <- sort(as.integer(sub("^sim", "", sim_dirs)))
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

# Returns "from--to" keys of detected and missed reticulations for one sim from detection_report.txt,
# the pairs that detected each, and the traced pairs; NULL when there is no report or no block for the sim.
read_detection <- function(sim_index) {
  if (is.null(REPORT_DIR)) return(NULL)
  report_file <- file.path(REPORT_DIR, "detection_report.txt")
  if (!file.exists(report_file)) return(NULL)
  lines <- readLines(report_file)
  start <- grep(paste0("^sim", sim_index, ":"), lines)
  if (length(start) != 1) return(NULL)
  rest <- lines[-seq_len(start)]
  n <- which(!grepl("^  ", rest))[1] - 1
  if (is.na(n)) n <- length(rest)

  detected <- character(0)
  missed <- character(0)
  detected_by <- character(0)
  traced <- list()
  for (line in rest[seq_len(n)]) {
    m <- regmatches(line, regexec("^  detected (\\S+)\\([a-z]+\\) -- (\\S+)\\([a-z]+\\)  by (.*)$", line))[[1]]
    if (length(m) == 4) {
      key <- paste(m[2], m[3], sep = "--")
      detected <- c(detected, key)
      pairs <- strsplit(m[4], ", ")[[1]]
      detected_by[key] <- paste(pairs, collapse = "\n")
      for (p in pairs) traced[[length(traced) + 1]] <- sub("\\([a-z]+\\)$", "", strsplit(p, "--")[[1]])
      next
    }
    m <- regmatches(line, regexec("^  missed   (\\S+) -- (\\S+)$", line))[[1]]
    if (length(m) == 3) missed <- c(missed, paste(m[2], m[3], sep = "--"))
  }
  keys <- vapply(traced, function(p) paste(sort(p), collapse = "--"), "")
  list(detected = detected, missed = missed, traced = traced[!duplicated(keys)], detected_by = detected_by)
}

# purple bracket right of the tip labels joining each traced tip pair, one column per pair
draw_traced_pairs <- function(traced, tip_names, shown_labels) {
  lp    <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  ntip  <- length(tip_names)
  x_end <- max(lp$xx[seq_len(ntip)]) + max(strwidth(shown_labels)) * 1.1
  step  <- diff(par("usr")[1:2]) * 0.03
  for (k in seq_along(traced)) {
    y  <- lp$yy[match(traced[[k]], tip_names)]
    x0 <- x_end + k * step
    segments(x_end, y, x0, y, col = "purple", lwd = 1.5)
    segments(x0, y[1], x0, y[2], col = "purple", lwd = 1.5)
  }
  legend("bottomright", legend = paste0("traced tip pair (", length(traced), ")"), col = "purple",
         lwd = 1.5, bty = "n", cex = 0.7)
}

# ── Read one sim's CSVs, plot it, and write the PDF ───────────────────────────

# kind = "" plots simN/simN_nodes.csv and simN_edges.csv (raw phy.G, export_csv);
# kind = "filtered_" plots simN/simN_filtered_nodes.csv and simN_filtered_edges.csv
# (the "no_hyb_nodes"-collapsed graph gene-tree enumeration actually runs on,
# export_filtered) -- same node IDs as the raw graph, just fewer nodes, so the
# two PDFs stay directly comparable.
plot_sim <- function(sim_index, kind = "") {
  nodes_file <- sim_csv_file(sim_index, kind, "nodes")
  edges_file <- sim_csv_file(sim_index, kind, "edges")

  nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
  edges <- read.csv(edges_file, stringsAsFactors = FALSE)

  phy <- build_phy(nodes, edges)
  internal_ids <- attr(phy, "internal_ids")

  # build_phy makes every node without tree children a tip, including lineages that went extinct but are
  # kept as a reticulation source (e.g. a donor whose continuation died). Mark those apart from living tips.
  tip_is_live <- as.character(nodes$is_leaf[match(phy$tip.label, nodes$label)]) %in% c("True", "TRUE", "true")
  tip_color   <- ifelse(tip_is_live, "black", "grey50")
  tip_names   <- phy$tip.label
  phy$tip.label <- ifelse(tip_is_live, phy$tip.label, paste(phy$tip.label, "†"))

  # ape::plot.phylo needs >= 2 tips, and a single tip has nothing worth showing, so skip it
  if (length(phy$tip.label) < 2) {
    cat("Skipping sim", sim_index, kind, "- only", length(phy$tip.label), "tip, nothing to plot\n")
    return(invisible(NULL))
  }

  out_file <- file.path(PLOT_DIR, paste0("sim", sim_index, "_", kind, "network.pdf"))
  pdf(out_file, encoding = "WinAnsi")  # WinAnsi has the dagger used to mark extinct lineages; the default doesn't
  on.exit(dev.off(), add = TRUE)

  title_label <- paste0("BDH Network - ", if (kind == "") "raw" else "filtered", " - sim", sim_index)

  # reticulations: red = detected in the report, grey dashed = missed,
  # blue = no report (same order as build_phy's phy$reticulation rows)
  ret_edges <- edges[edges$edge_type == "reticulation", ]
  detection <- if (kind == "filtered_") read_detection(sim_index) else NULL
  ret_col <- "blue"; ret_lty <- 1
  if (!is.null(detection) && nrow(ret_edges) > 0) {
    fwd <- paste(ret_edges$from, ret_edges$to, sep = "--")
    rev <- paste(ret_edges$to, ret_edges$from, sep = "--")
    hit  <- fwd %in% detection$detected | rev %in% detection$detected
    miss <- fwd %in% detection$missed   | rev %in% detection$missed
    ret_col <- ifelse(hit, "red", ifelse(miss, "grey50", "blue"))
    ret_lty <- ifelse(hit, 1, 2)
    title_label <- paste0(title_label, " - ", sum(hit), "/", nrow(ret_edges), " retic detected")
  }
  traced <- if (!is.null(detection)) Filter(function(p) all(p %in% tip_names), detection$traced) else list()
  depth  <- max(ape::node.depth.edgelength(phy))
  x_lim  <- if (length(traced) > 0) c(0, depth * (1.5 + 0.06 * length(traced))) else NULL
  plot(phy, main = title_label, col = ret_col, lty = ret_lty, alpha = 0.9, tip.color = tip_color,
       x.lim = x_lim)
  if (length(traced) > 0) draw_traced_pairs(traced, tip_names, phy$tip.label)
  if (!all(tip_is_live)) {
    legend("bottomleft", legend = "† not a living tip (extinct lineage)", text.col = "grey50", bty = "n", cex = 0.7)
  }
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
    # the pairs that detected the reticulation, under its length
    if (!is.null(detection) && length(detection$detected_by) > 0) {
      fwd <- paste(ret_edges$from, ret_edges$to, sep = "--")
      rev <- paste(ret_edges$to, ret_edges$from, sep = "--")
      leader_text <- ifelse(fwd %in% names(detection$detected_by), detection$detected_by[fwd],
                            ifelse(rev %in% names(detection$detected_by), detection$detected_by[rev], ""))
      text(mid_x, mid_y, leader_text, cex = 0.4, col = "red", adj = c(0.5, 1.4))
    }
  }
  if (!is.null(detection) && nrow(ret_edges) > 0) {
    legend("topleft", legend = c("retic detected", "retic missed"),
           col = c("red", "grey50"), lty = c(1, 2), bty = "n", cex = 0.7)
  }

  cat("Plot written to", out_file, "\n")
}

# ── Run ────────────────────────────────────────────────────────────────────────

for (sim_index in sim_indices) {
  for (kind in c("", "filtered_")) {
    # Not every pipeline exports both kinds, so skip silently when absent
    if (!file.exists(sim_csv_file(sim_index, kind, "nodes"))) next
    result <- tryCatch({
      plot_sim(sim_index, kind)
      TRUE
    }, error = function(e) {
      cat("Skipping sim", sim_index, "kind='", kind, "' - error:", conditionMessage(e), "\n")
      FALSE
    })
  }
}
