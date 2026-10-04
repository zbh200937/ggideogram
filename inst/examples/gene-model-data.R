# Resolve GFF3 transcript and feature relationships, retaining source attributes.
gene_model_data <- function(features) {
  transcripts <- features[features$Type == 'mRNA', , drop = FALSE]
  parents <- strsplit(ifelse(features$Type == 'mRNA', features$ID, features$Parent),
    ',', fixed = TRUE)
  models <- features[rep(seq_len(nrow(features)), lengths(parents)), , drop = FALSE]
  models$Transcript <- unlist(parents, use.names = FALSE)
  models <- models[models$Transcript %in% transcripts$ID, , drop = FALSE]
  parents <- strsplit(transcripts$Parent[match(models$Transcript, transcripts$ID)],
    ',', fixed = TRUE)
  models <- models[rep(seq_len(nrow(models)), lengths(parents)), , drop = FALSE]
  models$Gene <- unlist(parents, use.names = FALSE)
  rownames(models) <- NULL
  models
}
