

#' Shuffle, Partition, and Encode Datasets.
#'
#' Randomly shuffles the given dataset before seperating it into two sets which
#' are then encoded and ready for use in training a linear model.
#'
#' @param examples A list of aligned peptide sequences.
#'
#' @param labels The continuously distributed function for each of the given
#' examples.
#'
#' @param num_p1 The number of examples to put the in the first dataset. The
#' remainder will be put in the second dataset.
#'
#' @param encode Function to be used to encode the given examples. NULL by
#' default, returning the original unmodified peptide sequences.
#'
#' @export
shuffled_partitions <- function(examples, labels, num_p1, encode=NULL){
  num_examples <- length(labels)

  inds <- sample(1:num_examples, num_examples)

  examples <- examples[inds]
  labels <- labels[inds]

  p1_examples <- examples[1:num_p1]
  p1_labels <- labels[1:num_p1]
  p2_examples <- examples[(1 + num_p1):num_examples]
  p2_labels <- labels[(1 + num_p1):num_examples]

  if (! is.null(encode)){
    p1_examples <- encode(p1_examples)
    p2_examples <- encode(p2_examples)
  }

  rtrn <- list(p1_examples, p1_labels, p2_examples, p2_labels)
  names(rtrn) <- c('e1', 'l1', 'e2', 'l2')
  return(rtrn)
}

# Helper for encoding funcitons.
strs_to_matrix <- function(strs){
  assert_that(is.list(strs))

  strs_split <- lapply(strs, strsplit, split="")
  strs_split <- lapply(strs_split, unlist)
  strs_vecs <- as.vector(strs_split)

  strs_matrix <- do.call(rbind, strs_vecs)

  return(strs_matrix)
}


#' Encode Aligned Peptide Sequences into a One-hot Matrix
#'
#' @param AA_seqs List of aligned character sequences.
#'
#' @return Rank 3 array of peptide sequence encodings.
#'
#' @examples
#'
#' examples <- rhoData$data
#'
#' encode_onehot(examples)
#'
#' @export
#' @import assertthat
encode_onehot <- function(AA_seqs){
  num_AAs <- 20
  Amino_acids <- unlist(strsplit("ACDEFGHIKLMNPQRSTVWY", ''))
  AAs_to_num <- as.list(1:20)
  names(AAs_to_num) = Amino_acids

  assert_that(is.list(AA_seqs))

  AA_matrix <- strs_to_matrix(AA_seqs)

  num_seqs <- length(AA_seqs)
  num_residues <- nchar(AA_seqs[[1]])

  one_hot_array <- array(0L, dim= c(num_seqs, num_residues, num_AAs))

  for (x_ind in seq(num_seqs)){
    for (aa_ind in seq(num_residues)){
      AA <- toupper(AA_matrix[x_ind, aa_ind])
      if (AA == '-'){
        next
      }
      AA_num <- AAs_to_num[[AA]]
      one_hot_array[x_ind, aa_ind, AA_num] = 1
    }
  }

  return(one_hot_array)
}

#' Encode Aligned Peptide Sequences Using the VHSE8 Matrix
#'
#' @param AA_seqs List of aligned character sequences.
#'
#' @return Rank 3 array of peptide sequence encodings.
#'
#' @examples
#'
#' examples <- rhoData$data
#'
#' encode_physchem(examples)
#'
#' @export
#' @import Peptides
#' @import assertthat
encode_physchem <- function(AA_seqs){
  data("AAdata")
  scale <- AAdata$VHSE
  num_vhses <- 8

  get_VHSEs <- function(AA){
    res <- numeric(num_vhses)
    for (v_num in seq_along(res)){
      res[v_num] <- scale[[v_num]][[AA]]
    }
    return(res)
  }

  assert_that(is.list(AA_seqs))

  AA_matrix <- strs_to_matrix(AA_seqs)

  num_seqs <- length(AA_seqs)
  num_residues <- nchar(AA_seqs[[1]])

  vhse_array <- array(0L, dim= c(num_seqs, num_residues, num_vhses))

  for (x_ind in seq(num_seqs)){
    for (aa_ind in seq(num_residues)){
      AA <- toupper(AA_matrix[x_ind, aa_ind])
      if (AA == '-'){
        next
      }
      vhses <- get_VHSEs(AA)
      vhse_array[x_ind, aa_ind, 1:num_vhses] <- vhses
    }
  }

  return(vhse_array)
}

# [END]
