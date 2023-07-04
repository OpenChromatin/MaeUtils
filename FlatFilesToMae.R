# v0.0.0, in collaboration with Gpt4. Here, the comments were typed by OCG, then provided to Gpt4 in conjunction with a thorough description of the goal of the project. The code as it stands below is very likely non-functional.


"
In the identifyExperimentDataDfs function, I made an assumption about how the non-matching columns should be concatenated to set as row names. I used paste function to concatenate the non-matching columns element-wise and set them as row names. Please verify if this is the desired behavior.
In the CheckIfSummarizedExperiment function, I assumed that each element of the ExperimentList should be converted to a SummarizedExperiment object with the "counts" assay. Please confirm if this is the correct assumption.
In the loadInMetaData function, I assumed that each element of ListOfDfs represents a metadata object to be assigned to the MAE. Please verify if this is the intended behavior.
Please review these parts of the code and make any necessary adjustments based on your specific requirements.
"
FileNameToStem <- function(FilePaths) {
" This function takes a vector of absolute file paths, strips away the directory, strips away the file extension, then returns the remaining strings."
  fNames <- basename(FilePaths)
  fStems <- tools::file_path_sans_ext(fNames)
  return(fStems)
}

RawTabularFileLoader <- function(TabularFileName) {
"This function takes the filePath of a file that can be read into a data.frame, skips lines 
beginning with #, and assumes the files are assumed to have headers, then returns Df"
  return(tryCatch( read.table( TabularFileName, na.strings = c("NA", "-"), comment.char = "#", header = TRUE ), error = function(e) NULL ))
}

load_files_from_folder <- function(dataDir) {
" This function gets names of all files in a folder, removes directory names,  prepares a named list using names based on the files themsleves using the
 FileNameToStem() function. Next, the files are loaded as data.frames into the named list of data.frames using the function RawTabularFileLoader() ."
  subDirs <- list.dirs(dataDir)
  AllFiles <- list.files(dataDir, full.names = TRUE, include.dirs = FALSE)
  fPaths <- setdiff(subDirs[2:length(subDirs)], AllFiles)
  dfNames <- FileNameToStem(fPaths)
  dfList <- lapply(fPaths, RawTabularFileLoader)
  names(dfList) <- dfNames
  return(dfList)
}

findPrimaryIdByStringDistance <- function(data, coldata) {
" This function uses the stringdist::stringdistmatrix() function to calculate string distance between all candidate columns, and the header columns of the remaining files. 

This is done by first taking the headers of every df in the list, then making a table of the number of times all colnames in any header appears. 

Next, the first 20 columns of each data.frame in the list is read. Completely numeric columns are skipped. Second, if the length of the column read in is longer than the
length of the table of colnames, it is skipped, because the superset of all colnames should match all of the items of the true primaryId column at least once. To allow for
odd scenarios, though, columns up to 40% longer than the length of the table will be analyzed.

Columns not skipped will be added to a list called EligibleColumns, and then stringdist::stringdistmatrix()
will be used to compute best match for each element of the column successively from among all the strings in
the table (after which the element is paired with the colname, but removed from the tempTable, so it does not 
match more than one item in the column being considered). 

After pairing each item in the column EligibleColumns, the primaryId column is identifed as the column that had 
the smallest mean string distance between the matched pairs. 

As a final check to verify the accuracy of the selection, the exact difference between candidate pairs of 
strings will be identified using stringdiff(). If one or more identical changes can be made to all elements of 
the table of colnames and this results in an exact match, then the primaryId column is confirmed, and the 
data.frame containing that column is confirmed as the ColDataDf, which will be used to make the top level 
ColData object for the MAE. 

Initially, the findPrimaryIdByStringDistance function, the original code didn't include the logic to check for an exact match between the primaryId column and the colnames. I made an assumption and added the logic to assign an Inf match score if an exact match is found.

Finally, the function assigns primaryId to be the rownames of ColDataDf, and removes ColDataDf from the list of all data.frames. "

  match_scores <- numeric(length(data))
  for (i in 1:length(data)) {
    column <- data[[i]]
    
    if (any(coldata[[1]] %in% column)) {
      match_scores[i] <- Inf
    } else {
      dist_scores <- stringdist::stringdistmatrix(column, coldata[[1]])
      min_dist <- apply(dist_scores, 2, min)
      match_scores[i] <- sum(1 / min_dist)
    }
  }
  ColDataDf <- coldata[[which.max(match_scores)]]
  colnames(ColDataDf) <- NULL
  return(list(ColDataDf = ColDataDf))
}

assignDfColClasses <- function(dfList) {
"Every column of each data.frame whose header name matches an element of primaryId will be assigned a class as follows. First, 
entirely numeric columns (after temporary removal of NAs) will be  assigned a Class as.integer() if 
all( is.integer(colElements) ) == TRUE, or as.numeric() otherwise. Columns containing either all strings or a mixture of strings 
and numerics will be assigned to characters. 

If the column is already a factor on import, move on; in either other event (numeric or character), next the RepeatedValuesToFactor() 
function is called. After calling unique() on the col, the quotient of the number of unique elements / the total length of the column 
is computed. If this quotient is < 0.25, then the column is re-assigned to a factor. The reference level of the factor is designated 
as either the lowest number (if the col was orginally designated numeric) or the most common string (if originally designated a character).

The function then returns the column either in its original format or as a factor."
  df <- as.data.frame(lapply(dfList, function(col_data) {
    if (all(sapply(col_data, is.numeric))) {
      if (all(sapply(col_data, function(x) is.integer(x) && !is.na(x)))) {
        return(as.integer(col_data))
      } else {https://github.com/OpenChromatin/MaeUtils/tree/main
        return(as.numeric(col_data))
      }
    } else if (all(sapply(col_data, is.character))) {
      col_data <- as.character(col_data)
      table_ratio <- length(table(col_data)) / length(col_data)
      if (table_ratio < 0.2) {
        if (is.numeric(col_data)) {
          min_value <- min(col_data, na.rm = TRUE)
          col_data <- as.factor(col_data)
          levels(col_data) <- c(as.character(min_value), levels(col_data)[-1])
        } else {
          freq_table <- table(col_data)
          most_freq_string <- names(freq_table)[which.max(freq_table)]
          col_data <- as.factor(col_data)
          levels(col_data) <- c(most_freq_string, levels(col_data)[-1])
        }
      }
      return(col_data)
    }
    return(col_data)
  }))
  return(df)
}

identifyExperimentDataDfs <- function(ListOfDfs, ColDataDf) {
"This function compares the primaryId column to the colnames of each data Df in turn. Any element in colnames that does not 
exactly match an element of the primaryId is appended to a vector, then these columns are pasted together element-wise and 
set as the row.names() for each Df.

The remaining columns (that have exact matches) will remain in the data object. Finally, the 
function will remove the Dfs identified as Experiment Dfs from ListOfDfs, and return both." 
  remaining_cols <- character()
  for (i in seq_along(ListOfDfs)) {
    cols <- colnames(ListOfDfs[[i]])
    non_matching_cols <- setdiff(cols, ColDataDf)
    remaining_cols <- c(remaining_cols, non_matching_cols)
    row_names <- apply(ListOfDfs[[i]][remaining_cols], 1, paste, collapse = "_")
    row.names(ListOfDfs[[i]]) <- row_names
    ListOfDfs[[i]] <- ListOfDfs[[i]][, ColDataDf]
  }
  return(list(ExperimentList = ListOfDfs, RemainingDfs = remaining_cols))
}

CheckIfSummarizedExperiment <- function(ExperimentList) {
"This fnuction checks whether each Df in ExperimentList can be converted to a summarized experiment. 
If so, this will be undertaken, and the ColDataDf will be subsetted to just the fraction pertaining
to the current Df. This will then be set as colData of the summarizedExperiment, which will ultimately
be rolled into the MAE. 

Finally, the function then returns the object, either unaltered or as a SE, back to a (now updated) ExperimentList."
  UpdatedExperimentList <- list()
  for (i in seq_along(ExperimentList)) {
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = ExperimentList[[i]]))
    UpdatedExperimentList[[i]] <- se
  }
  return(UpdatedExperimentList)
}

loadInMetaData <- function(ListOfDfs) {
"Any Df that remains in list of Dfs corresponds to MetaData, and can be assigned to the MAE object using metadata(MAE)<-object."
  MetaDataObjects <- list()
  for (i in seq_along(ListOfDfs)) { MetaDataObjects[[i]] <- ListOfDfs[[i]] }
  return(MetaDataObjects)
}

FlatFilesToMae <- function(DataDir) {
  AllInputDfs <- load_files_from_folder(DataDir)
  AllInputDfs <- assignDfColClasses(AllInputDfs)
  results <- findPrimaryIdByStringDistance(AllInputDfs, AllInputDfs[[1]])
  ColDataDf <- results$ColDataDf
  RemainingInputDfs <- AllInputDfs[-1]
  ExpResults <- identifyExperimentDataDfs(RemainingInputDfs, ColDataDf)
  ExperimentList <- ExpResults$ExperimentList
  RemainingDfs <- ExpResults$RemainingDfs
  ExperimentList <- CheckIfSummarizedExperiment(ExperimentList)
  MetaDataToAttach <- loadInMetaData(RemainingDfs)
  
  MAE <- MultiAssayExperiment::MultiAssayExperiment(
    colData = ColDataDf,
    experiments = ExperimentList,
    metadata = MetaDataToAttach
  )
  
  return(MAE)
}


FlatFilesToMae<-function(DataDir) {
"This is the master function; it calls
	LoadFilesInFolder()
	assignDfColClasses()
	findColDataAndPrimaryId()
	identifyExperimentDataDfs()
	checkIfSummarizedExperiment
	identifyMetaData()
	sequentially; the remaining functions are called from within these. 

	Finally, the function returns the MAE."
#
	libsNeeded<-c('stringdist','SummarizedExperiment','MultiAssayExperiment')
	BiocManager::install(libsNeeded, update = TRUE, ask = FALSE)
	lapply(libsNeeded, library, character.only = TRUE)
#
	AllInputDfs<-LoadFilesInFolder(DataDir)
	AllInputDfs<-assignDfColClasses(DataDir)
	results<-findColDataAndPrimaryId(DataDir); ColDataDf<-results$ColDataDf; RemainingInputDfs<-results$RemainingInputDfs
	ExpResults<-makeExperimentList(RemainingInputDfs, ColDataDf); ExperimentList<-ExpResults$ExperimentList; RemainingDfs<-ExpResults$RemainingDfs
	ExperimentList<-checkIfSummarizedExperiment(ExperimentList)
	MetaDataToAttach<-identifyMetaData(RemainingDfs)
#
	MAE<-MultiAssayExperiment(colData=ColDataDf, experiments=ExperimentList, metadata=MetaDataToAttach)
	return(MAE)
}

MultiAssayExperiment(DirWithFiles)
