.libPaths("~/R/library")
library(ape)
library(SiPhyNetwork)

# ── Parameters ────────────────────────────────────────────────────────────────

OUT_DIR   <- "../../Outputs/phylo_outputs"
SIM_INDEX <- 1   # change to plot a different simulation

nodes <- read.csv(file.path(OUT_DIR, paste0("sim", SIM_INDEX, "_nodes.csv")),
                  stringsAsFactors = FALSE)
edges <- read.csv(file.path(OUT_DIR, paste0("sim", SIM_INDEX, "_edges.csv")),
                  stringsAsFactors = FALSE)

# ── Build ape-compatible node IDs ─────────────────────────────────────────────

tree_edges <- edges[edges$edge_type == "tree", ]
ret_edges  <- edges[edges$edge_type == "reticulation", ]

# Tips: nodes with no outgoing tree edges
tip_ids      <- nodes$id[!(nodes$id %in% tree_edges$from)]
# Root: node with no incoming edges (not a child in any edge)
root_id      <- nodes$id[!(nodes$id %in% c(tree_edges$to, ret_edges$to))]
# Internal: everything else (put root first so it maps to ntips + 1)
internal_ids <- c(root_id, setdiff(nodes$id[nodes$id %in% tree_edges$from], root_id))

ntips <- length(tip_ids)
Nnode <- length(internal_ids)

# Mapping: sim node ID -> ape node ID
id_map <- c(setNames(seq_len(ntips),           as.character(tip_ids)),
            setNames(seq(ntips + 1, ntips + Nnode), as.character(internal_ids)))

# ── Construct phylo object ────────────────────────────────────────────────────

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

# ── Add reticulation edges ────────────────────────────────────────────────────

if (nrow(ret_edges) > 0) {
  phy$reticulation <- matrix(
    c(id_map[as.character(ret_edges$from)],
      id_map[as.character(ret_edges$to)]),
    ncol = 2)
  phy$reticulation.length <- ret_edges$length
  class(phy) <- c("evonet", "phylo")
}

# ── Plot ──────────────────────────────────────────────────────────────────────

out_file <- file.path(OUT_DIR, paste0("sim", SIM_INDEX, "_network.pdf"))
pdf(out_file)
plot(phy, main = paste("BDH Network - sim", SIM_INDEX))

internal_labels <- nodes$label[match(internal_ids, nodes$id)]
ape::nodelabels(text = internal_labels, cex = 0.6, frame = "none", col = "blue")

dev.off()
cat("Plot written to", out_file, "\n")
