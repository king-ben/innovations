setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(XML)
library(devtools)
source_url("https://raw.githubusercontent.com/king-ben/beast_asr/main/asr_functions.R")

wp <- find_partitions("Malekula.nex")
wn <- find_site_names("Malekula.nex")
links <- rate_links(wp, cutoff=5)

#innovations_blocks function
# wp - the positions for each partition
# wn - the names of each cognate in the order they appear in the dataset
# links - the name of the sitemodel for each partition when using binning
# taxonset - the name of a taxonset that should be established somewhere in the xml file. innovations to the common ancestor of this set are logged
# model - preset to avoid having to manually input the from and to states, at the moment either covarion or standard
# fromstates and tostates - can manually put in any states to log between if the presets are not desired. useful e.g. for logging losses instead of gains
# clockname - the id for the clock if not RelaxedClock.c:clock
# logevery - how often to log. default 1000

dist <- innovation_blocks(wp, wn, links, taxonset="east", model="covarion", clockname="@OptimisedRelaxedClock.c:clock", logevery=10000)
saveXML(dist, file="innovation_blocks.xml")