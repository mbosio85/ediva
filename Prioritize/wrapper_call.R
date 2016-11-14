thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}


args <- commandArgs(TRUE)
model <- as.character(args[1])
file_location <- as.character(args[2])
outfile <- as.character(args[3])
script.dir <-dirname(thisFile())
source(paste(script.dir,'/predict_model_file.R',sep=''))

predict_me(model,file_location,outfile)


