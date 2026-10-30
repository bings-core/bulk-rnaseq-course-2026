# Build the TxDb that 02_deg.R reads (gene and transcript models) from the GENCODE GTF.
# Run by setup_resources.sh inside the R container; the file name matches the one 02_deg.R loads.
suppressPackageStartupMessages(library(GenomicFeatures))

species <- Sys.getenv("SPECIES")
ref     <- file.path(Sys.getenv("RESOURCE_DIR"), species)
name    <- if (species == "human") "homo_sapiens_grch38_gencode_36" else "mus_musculus_grcm38_gencode_M25"

txdb <- makeTxDbFromGFF(file.path(ref, "gtf", "annotation.gtf"), format = "gtf")
saveDb(txdb, file.path(ref, "txdb", paste0(name, ".sqlite")))
cat("TxDb written to", file.path(ref, "txdb", paste0(name, ".sqlite")), "\n")
