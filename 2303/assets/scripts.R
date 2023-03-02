# Additional R scripts

#' @title download_data
#' @description Used for externally imported objects. Checks if the object is available, else downloads it from the github master repo.
#' @param fileName Name of the object to download. Example: count_raw.txt
#' @param directory Name of the directory where the object should be downloaded. Example: results/06-R-analyses 
#' @param dirRepo Directory where object is found in the repo. Example: data/ (default)
#' @param repo A full github repo name. Example: NBISweden/workshop-RNAseq (default)
#' @param branch Github repo branch name. Defaults to "master".
#' @example
#' 
#' download_data(fileName = "metadata_raw.txt", directory = "results/06-R-analyses", repo = "NBISweden/workshop-RNAseq")
#' 
#' 
download_data2 <- function(fileName = NULL,
                          directory = NULL, 
                          dirRepo = "data", 
                          repo = "NBISweden/workshop-RNAseq", 
                          branch = "master"){

  if(!is.null(directory) && !is.null(repo)){
    directory <- gsub("/$", "", gsub("^/|^./", "", directory))
    dirRepo <- gsub("/$", "", gsub("^/|^./", "", dirRepo))
    
    # if file doesn't exist
    if(!file.exists(paste0(directory, "/", fileName))) {
      
      # if directory doesn't exist, create it
      if(!dir.exists(directory)) {
        message(paste0("Creating directory: ", directory))
        dir.create(directory, recursive = TRUE)
      }
      
      download.file(url = 
                      paste0("https://raw.githubusercontent.com/",
                             repo,"/",
                             branch,"/",
                             paste0(dirRepo, "/", fileName)), 
                    paste0(directory, "/", fileName))
    }
  }
}

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