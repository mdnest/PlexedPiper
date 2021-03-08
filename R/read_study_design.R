#' Utilities for reading study design files.
#'
#' @description
#' Fetches study design results from either local folder of DMS.
#' Checks that study design files are internally consistent.
#' DMS functionality not useful outside of PNNL unless connected through VPN.
#'
#' @description
#' * `read_study_design()`: returns a list of study design tables, accessible by $
#'
#' @param path_to_folder (string) path to folder containing study design files
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr filter
#'
#' @name read_study_design
#'
#' @examples
#' study_design <- read_study_design("study_design_folder')
#' fractions <- study_design$fractions
#' samples <- study_design$samples
#' references <- study_design$references


#' @export
#' @rdname read_study_design
# gets 3 study design files from local directory
read_study_design <- function(path_to_study_design) {
  
  
  ## fetch samples.txt
  pathToFile <- list.files(path=path_to_study_design,
                           pattern="^samples.txt$",
                           full.names=T)
  if(length(pathToFile) == 0) {
    stop("'samples.txt' not found.")
  }
  
  samples <- readr::read_tsv(pathToFile,
                             col_types=readr::cols(.default = "c"),
                             progress=FALSE)
  if (!setequal(colnames(samples), c("PlexID",
                                     "QuantBlock",
                                     "ReporterAlias",
                                     "ReporterName",
                                     "MeasurementName"))) {
    stop("There are incorrect column names or missing columns in the 'samples'
         study design table.")
  }
  
  ## fetch fractions.txt
  pathToFile <- list.files(path=path_to_study_design,
                           pattern="^fractions.txt$",
                           full.names=T)
  if (length(pathToFile) == 0){
    stop("'fractions.txt' not found.")
  }
  
  
  fractions <- readr::read_tsv(pathToFile,
                               col_types=readr::cols(.default = "c"),
                               progress=FALSE)
  if (!setequal(colnames(fractions), c("PlexID",
                                       "Dataset"))) {
    stop("There are incorrect column names or missing columns in the 'fractions'
         study design table.")
  }
  
  ## fetch references.txt
  pathToFile <- list.files(path=path_to_study_design,
                           pattern="^references.txt$",
                           full.names=T)
  if(length(pathToFile) > 0) {
    references <- readr::read_tsv(pathToFile,
                                  col_types=readr::cols(.default = "c"),
                                  progress=FALSE)
    if (!setequal(colnames(references), c("PlexID",
                                          "QuantBlock",
                                          "Reference"))) {
      stop("There are incorrect column names or missing columns in the 'references'
           study design table.")
    }
  } else {
    warning("'references.txt' not found. It will be made automatically from `samples.txt`")
    references <- filter(samples, is.na(MeasurementName)) %>%
      dplyr::select(PlexID, QuantBlock, ReporterAlias) %>%
      dplyr::rename(Reference = ReporterAlias)
  }
  
  # Check for duplicates
  if (any(duplicated(fractions$Dataset))) {
    stop("Duplicate datasets in 'fractions.txt'")
  }
  if (any(duplicated(samples$MeasurementName[!is.na(samples$MeasurementName)]))) {
    stop("Duplicate sample names in 'samples.txt'")
  }
  if (!setequal(fractions$PlexID, samples$PlexID)) {
    stop("Plex IDs in 'fractions.txt' and 'samples.txt' do not match.")
  }
  # samples_subset_to_references <- filter(samples, is.na(MeasurementName)) %>%
  #   dplyr::select(PlexID, QuantBlock, ReporterAlias) %>%
  #   dplyr::rename(Reference = ReporterAlias)
  # if (!identical(samples_subset_to_references, references)) {
  #   stop("'samples.txt' and 'references.txt' are incompatible.")
  # }
  
  study_design <- list(samples = samples,
                    fractions = fractions,
                    references = references)
  
  return(study_design)
}