#' Format Tables for BIC
#'
#' @description
#' Assembles data in format compliant with BIC requirements. 
#'
#' @description
#' * `make_rii_peptide_gl()`: returns 'RII_peptide.txt' table formatted for BIC (global)
#' * `make_results_ratio_gl()`: returns 'results_ratio.txt' table (global)
#' * `make_rii_peptide_ph()`: returns 'RII_peptide.txt' table (phospho)
#' * `make_results_ratio_ph()`: returns 'results_ratio.txt' table (phospho)
#' * `map_site_flanks()`: returns MSnID object with +/- 7 amino acids (by default) surrounding each PTM
#' @md
#'
#' @param msnid (MSnID-object) final filtered version of MSnID object
#' @param masic_data (data.frame) final filtered version of MASIC table
#' @param fractions (data.frame) Study design table linking Dataset with PlexID
#' @param samples (data.frame) Study design table linking sample names with TMT channels and PlexID
#' @param references (data.frame) Study design table describing reference value calculation
#' @param org_name (character) Organism name. Default is 'Rattus norvegicus'
#' @param sep (character) Single character used to concatenate protein, SiteID, and peptide
#'
#' @importFrom dplyr select inner_join mutate %>% case_when rename group_by summarize
#' @importFrom tibble rownames_to_column
#' @importFrom purrr pmap
#'
#'
#' @name motrpac_bic_output
#'
#' @examples
#' # Prepare MS/MS IDs
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
#' show(msnid)
#' path_to_FASTA <- system.file("extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", package = "PlexedPiperTestData")
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' msnid <- filter_msgf_data_protein_level(msnid, 0.01)
#' show(msnid)
#'
#' # Prepare table with reporter ion intensities
#' path_to_MASIC_results <- system.file("extdata/global/masic_output", package = "PlexedPiperTestData")
#' masic_data <- read_masic_data(path_to_MASIC_results, extra_metrics=TRUE)
#' masic_data <- filter_masic_data(masic_data, 0.5, 0)
#' 
#' library(readr)
#' fractions <- read_tsv(system.file("extdata/study_design/fractions.txt", package = "PlexedPiperTestData"))
#' samples <- read_tsv(system.file("extdata/study_design/samples.txt", package = "PlexedPiperTestData"))
#' references <- read_tsv(system.file("extdata/study_design/references.txt", package = "PlexedPiperTestData"))
#' 
#' results_ratio <- make_results_ratio_gl(msnid, masic_data, fractions, samples, references, org_name = "Rattus norvegicus")
#' rii_peptide <- make_rii_peptide_gl(msnid, masic_data, fractions, samples, references, org_name = "Rattus norvegicus")


#' @export
#' @rdname motrpac_bic_output
make_rii_peptide_gl <- function(msnid, masic_data, fractions, samples, 
                                references, org_name = "Rattus norvegicus") {
  samples_rii <- samples %>%
    mutate(MeasurementName = case_when(is.na(MeasurementName) ~ paste0("Ref", "_", PlexID),
                                       TRUE ~ MeasurementName))
  
  
  references_rii <- references %>%
    mutate(Reference = 1)
  
  aggregation_level <- c("accession", "peptide")
  crosstab <- create_crosstab(msnid, 
                              masic_data, 
                              aggregation_level, 
                              fractions, samples_rii, references_rii)
  
  crosstab <- 2^crosstab # undo log2
  
  x <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  y <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID")
  conv <- inner_join(x, y) %>% 
    dplyr::rename(protein_id = REFSEQ,
                  gene_symbol = SYMBOL,
                  entrez_id = ENTREZID)
  
  crosstab <- crosstab %>% 
    as.data.frame() %>% 
    rownames_to_column('ids') %>% 
    mutate(protein_id = sub("(^.*)\\.\\d+@.*", "\\1", ids),
           sequence = sub("(^.*)@(.*)", "\\2", ids),
           organism_name = org_name) %>% 
    dplyr::select(-ids) %>%
    inner_join(conv) %>% 
    dplyr::select(protein_id, organism_name, sequence, gene_symbol, entrez_id, everything())
  
  return(crosstab)
  
}


#' @export
#' @rdname motrpac_bic_output
make_results_ratio_gl <- function(msnid, masic_data, fractions, samples, 
                                  references, org_name = "Rattus norvegicus") {
  aggregation_level <- c("accession")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)
  
  x <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  y <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID")
  conv <- inner_join(x, y) %>% 
    dplyr::rename(protein_id = REFSEQ,
                  gene_symbol = SYMBOL,
                  entrez_id = ENTREZID)
  
  crosstab <- crosstab %>% 
    as.data.frame() %>% 
    rownames_to_column('protein_id')
  
  crosstab <- psms(msnid) %>%
    dplyr::select(accession, peptide) %>%
    group_by(accession) %>%
    summarize(num_peptides = n()) %>%
    dplyr::rename(protein_id = accession) %>%
    left_join(crosstab, .)
  
  crosstab <- crosstab %>%
    dplyr::mutate(protein_id = sub("(.P_.*)\\.\\d+", "\\1", protein_id),
                  organism_name = org_name)
  
  results_ratio <- inner_join(crosstab, conv) %>% 
    dplyr::select(protein_id, organism_name, gene_symbol, entrez_id, num_peptides,
                  everything())
  
  return(results_ratio)
}

#' @export
#' @rdname motrpac_bic_output
make_rii_peptide_ph <- function(msnid, masic_data, fractions, samples, references,
                                org_name = "Rattus norvegicus", sep="_") {
  samples_rii <- samples %>%
    mutate(MeasurementName = case_when(is.na(MeasurementName) ~ paste0("Ref", "_", PlexID),
                                       TRUE ~ MeasurementName))
  
  
  references_rii <- references %>%
    mutate(Reference = 1)
  
  aggregation_level <- c("accession", "peptide")
  crosstab <- create_crosstab(msnid, 
                              masic_data, 
                              aggregation_level, 
                              fractions, samples_rii, references_rii)
  
  crosstab <- 2^crosstab # undo log2
  crosstab <- crosstab %>% 
    as.data.frame() %>% 
    rownames_to_column("Specie") %>% 
    mutate(protein_id = sub("(^.*)\\.\\d+@.*", "\\1", Specie),
           sequence = sub("(^.*)@(.*)", "\\2", Specie)) %>% 
    dplyr::select(protein_id, sequence, everything(), -Specie)
  
  conv <- dplyr::select(psms(msnid),
                        ptm_id = SiteID,
                        protein_id = accession,
                        sequence = peptide,
                        VMsiteFlanks) %>% 
    mutate(protein_id = sub("(.P_.*)\\.\\d+$", "\\1", protein_id))
  
  crosstab <- left_join(crosstab, conv) %>% 
    distinct() %>% 
    mutate(ptm_peptide = paste0(ptm_id, "-", sequence))
  
  ## Add Genes + EntrezID
  x <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  y <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID")
  conv <- inner_join(x, y) %>% 
    dplyr::rename(protein_id = REFSEQ,
                  gene_symbol = SYMBOL,
                  entrez_id = ENTREZID)
  crosstab <- inner_join(crosstab, conv)
  ## Add A-score
  ascore <- dplyr::select(psms(msnid), protein_id = accession, 
                          sequence = peptide, confident_score = maxAScore) %>% 
    mutate(protein_id = sub("(.P_.*)\\.\\d+", "\\1", protein_id)) %>% 
    group_by(protein_id, sequence) %>% 
    summarize(confident_score = max(confident_score))
  
  crosstab <- left_join(crosstab, ascore) %>% 
    mutate(confident_site = dplyr::case_when(confident_score >= 17 ~ TRUE,
                                             confident_score < 17 ~ FALSE),
           ptm_id = gsub("-", sep, ptm_id),
           ptm_peptide = gsub("-", sep, ptm_peptide),
           organism_name = org_name) %>% 
    dplyr::select(protein_id, organism_name, sequence, ptm_id, ptm_peptide, VMsiteFlanks, gene_symbol, entrez_id, confident_score, confident_site, everything()) %>% 
    distinct()
  return(crosstab)
}

#' @export
#' @rdname motrpac_bic_output
make_results_ratio_ph <- function(msnid, masic_data, fractions, samples, 
                                  references, org_name = "Rattus norvegicus", sep="_") {
  
  aggregation_level <- c("SiteID")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)
  crosstab <- crosstab %>% 
    as.data.frame() %>% 
    rownames_to_column('ptm_id') %>% 
    mutate(protein_id = sub("(.P_.*)\\.\\d+-.*", "\\1", ptm_id))
  
  x <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  y <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID")
  conv <- inner_join(x, y) %>% 
    dplyr::rename(protein_id = REFSEQ,
                  gene_symbol = SYMBOL,
                  entrez_id = ENTREZID)
  crosstab <- inner_join(crosstab, conv)
  
  ascore <- dplyr::select(psms(msnid), protein_id = accession, 
                          ptm_id = SiteID, 
                          VMsiteFlanks,
                          confident_score = maxAScore) %>% 
    mutate(protein_id = sub("(.P_.*)\\.\\d+", "\\1", protein_id)) %>% 
    group_by(protein_id, ptm_id, VMsiteFlanks) %>% 
    summarize(confident_score = max(confident_score))
  
  crosstab <- left_join(crosstab, ascore)
  
  crosstab <- crosstab %>% 
    mutate(confident_site = dplyr::case_when(confident_score >= 17 ~ TRUE,
                                             confident_score < 17 ~ FALSE),
           ptm_id = gsub("-", sep, ptm_id),
           organism_name = org_name) %>% 
    dplyr::select(ptm_id, VMsiteFlanks, organism_name, protein_id, gene_symbol, entrez_id, confident_score, confident_site, everything())
  
  crosstab[, c(9:ncol(crosstab))] <- signif(crosstab[, c(9:ncol(crosstab))], 3)
  return(crosstab)
  
}












#' @export
#' @rdname motrpac_bic_output
map_site_flanks <- function (msnid, radius=7L, collapse="|") {
  
  x <- psms(msnid) %>%
    dplyr::select(Peptide, pepSeq, ModAAs, ModShift) %>%
    distinct()
  
  f <- function(pepSeq_i, ModAA_i, ModShift_i) {
    VMsiteFlanks <- c()
    for (k in unlist(ModShift_i)) {
      k1 <- k+1 # possible bug. why is k off by 1?
      site_left <- substr(pepSeq_i, max(k1-radius,1), k1-1)
      site_right <- substr(pepSeq_i, k1+1, min(k1+radius,nchar(pepSeq_i)))
      
      mod_aa <- tolower(substr(pepSeq_i, k1, k1))
      
      VMsite <- paste0(site_left, tolower(mod_aa), site_right)
      VMsiteFlanks <- c(VMsiteFlanks, VMsite)
    }
    VMsiteFlanks <- paste(VMsiteFlanks, collapse=collapse)
    return(VMsiteFlanks)
  }
  
  
  x$VMsiteFlanks <-  pmap(list(x$pepSeq, x$ModAAs, x$ModShift), f)     
  
  msnid@psms <- psms(msnid) %>% mutate(VMsiteFlanks=NULL) %>%
    left_join(x) %>% data.table()
  
  return(msnid)
}



#' @export
#' @rdname motrpac_bic_output
assess_redundant_accessions <- function(msnid, collapse="|") {
  msnid@psms <- psms(msnid) %>%
    dplyr::select(accession, peptide) %>%
    distinct() %>%
    group_by(peptide) %>%
    summarize(reduntant_accessions = paste(accession, collapse=collapse)) %>%
    left_join(psms(msnid), .) %>%
    data.table()
  return(msnid)
}


