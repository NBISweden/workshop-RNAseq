# Additional R scripts

#' @title download_data
#' @description Used for externally imported objects. Checks if the object is available, else downloads it from the github master repo.
#' @param path Full path to object. Example: data/count_raw.txt
#' @param repo A full github repo name. Example: NBISweden/workshop-RNAseq
#' @param branch Github repo branch name. Defaults to "master".
#' @example
#' 
#' download_data("data/metadata_raw.txt","NBISweden/workshop-RNAseq")
#' 
download_data <- function(path=NULL,repo="NBISweden/workshop-RNAseq",branch="master"){

  if(!is.null(path) && !is.null(repo)){
    path <- gsub("^/|^./","",path)
    
    # if path doesn't exist
    if(!file.exists(path)) {
      
      # if directory doesn't exist, create it
      if(!dir.exists(dirname(path))) {
        message(paste0("Creating directory: ",dirname(path)))
        dir.create(dirname(path))
      }
      
      download.file(url=paste0("https://raw.githubusercontent.com/",repo,"/",branch,"/",path),path)
    }
  }
}


