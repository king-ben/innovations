setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(data.table)
log <- fread("asr_logger.txt")
log <- log[round(0.1*nrow(log)):nrow(log),]
probs <- apply(log, 2, mean)
probs <- probs[order(probs, decreasing =TRUE)]
innovations <- probs[probs>0.5]
length(innovations)