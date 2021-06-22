library(rotl)
library(phylocomr)
# devtools::install_github("fmichonneau/phyloch")
library(phyloch)
# devtools::install_github("phylotastic/datelife")
library(datelife)
library(phytools)
data(strat2012)

## Load file with correspondence between OTL and pollen taxonomy
taxa <- read.csv2("Data/Pollen_taxa_changes.csv")

# Query data to the OTL ddbb
resolved_names <- tnrs_match_names(as.character(taxa$pollen_otl_names), context_name = "Seed plants")
                                   
# Download tree topology
pollen_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id, label_format = "name")

# Download age data from datelife project
pollen_tree_dr <- get_datelife_result(pollen_tree)

# Summarzie chronograms
pollen_tree_phylo <-  summarize_datelife_result(pollen_tree_dr, summary_format = "phylo_median")

# Calibrate tree based in existing chronograms
pollen_tree_calibs <- get_all_calibrations(pollen_tree_phylo)

# Use calibrations to estimate ages for not calibrated nodes
pollen_tree_bladj <- use_calibrations(pollen_tree, pollen_tree_calibs, "bladj")

# Remove names suffixes included by OTL for some taxa
new.tip.labels <- strsplit(pollen_tree_bladj$tip.label, "_\\-") # In Linux OS  
# new.tip.labels <- strsplit(pollen_tree_bladj$tip.label, "_\\(") # In Windows OS 
new.tip.labels <- lapply(new.tip.labels, function(x){x[1]})
pollen_tree_bladj$tip.label <- unlist(new.tip.labels)
rm(new.tip.labels)

# Change names by the corresponding pollen taxa names
order <- match(pollen_tree_bladj$tip.label, taxa$otl_names)
pollen_tree_bladj$tip.label <- as.character(taxa$pollen_names[order])

# Remove gymnosperms from the tree
gymnosperm_genus <- c("Juniperus.Thuja", "Taxodium", "Taxus", "Ephedra", "Pinus", "Picea", "Larix", "Abies", "Tsuga")
gymnosperm_pollen_tree <- keep.tip(pollen_tree_bladj, gymnosperm_genus)
angiosperm_pollen_tree <- drop.tip(pollen_tree_bladj, gymnosperm_genus)

# Save file in the hard-disk drive
write.tree(pollen_tree_bladj, "Results/Data/pollen_tree.nwk")
write.tree(gymnosperm_pollen_tree, "Results/Data/gymnosperm_pollen_tree.nwk")
write.tree(angiosperm_pollen_tree, "Results/Data/angiosperm_pollen_tree.nwk")


# Plot phylogenetic tree with node ages and geologic scale
phy <- pollen_tree_bladj
phy$tip.label[c(8, 55, 59, 68, 72)] <- c("Planera", "Ambrosia", "Linnaea", "Chamaedaphne","Cyrilla")

# Prepare ages from callibration nodees for plotting
nodes.ages <- round(phy$used_calibrations, 0)
nodes.names <- names(nodes.ages)
nodes.numbers <- which(phy$node.label %in% nodes.names)

# Ploting tree
pdf(file="Results/Figures/OToL.pdf", width=5, height=8)
  plot(phy,
       y.lim=c(-30, 110),
       cex=0.35,
       edge.color=c(rep("blue", 192), rep("darkgreen", 17)))
  nodelabels(text=nodes.ages,
             node=nodes.numbers + 106,
             cex=0.4,
             )
  geo.legend()
dev.off()
