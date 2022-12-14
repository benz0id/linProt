#' Shuffle, Partition, and Encode Datasets.
#'
#' Randomly shuffles the given dataset before seperating it into two sets which
#' are then encoded and ready for use to train a linear model.
#'
#' @param examples A list of aligned amino acid sequences.
#'
#' @param labels A list of floating point numbers representing the
#' experimentally derived protein function for each example sequence.
#'
#' @param num_p1 The number of examples to put the in the first dataset. The
#' remainder will be put in the second dataset. This should be greater than 0
#' and less than or equal to the number of eaxamples in the given dataset.
#'
#' @param encode Function to be used to encode the given examples. NULL by
#' default, returning the original unmodified peptide sequences. The two
#' encoding functions provided by this package are `encode_onehot` and
#' `encode_physchem`.
#'
#'
#' @examples
#' examples <- rhoData$data
#' labels <- rhoData$labels
#'
#' # Shuffle, partition, and encode data set.
#' shuffled_datasets <- shuffled_partitions(examples, labels, 650,
#'                                          encode=encode_onehot)
#'
#' train_data <- shuffled_datasets$e1
#' train_labels <- shuffled_datasets$l1
#'
#' valid_data <- shuffled_datasets$e2
#' valid_labels <- shuffled_datasets$l2
#'
#'
#' @export
shuffled_partitions <- function(examples, labels, num_p1, encode=NULL){
  num_examples <- length(labels)

  # Shuffle the data.
  inds <- sample(1:num_examples, num_examples)

  examples <- examples[inds]
  labels <- labels[inds]

  # Partition the data.
  p1_examples <- examples[1:num_p1]
  p1_labels <- labels[1:num_p1]
  p2_examples <- examples[(1 + num_p1):num_examples]
  p2_labels <- labels[(1 + num_p1):num_examples]

  # Encode the data.
  if (! is.null(encode)){
    p1_examples <- encode(p1_examples)
    p2_examples <- encode(p2_examples)
  }

  rtrn <- list(p1_examples, p1_labels, p2_examples, p2_labels)
  names(rtrn) <- c('e1', 'l1', 'e2', 'l2')
  return(rtrn)
}

# Helper for encoding functions.
strs_to_matrix <- function(strs){
  assert_that(is.list(strs))
  # Convert list of strings to matrix of characters.
  strs_split <- lapply(strs, strsplit, split="")
  strs_split <- lapply(strs_split, unlist)
  strs_vecs <- as.vector(strs_split)

  strs_matrix <- do.call(rbind, strs_vecs)

  return(strs_matrix)
}


#' Encode a list of aligned AA sequences into a tensor of one-hot encodings.
#'
#' @param AA_seqs A list of aligned peptides. Represented as a list of character
#' vectors.
#'
#' @return Rank 3 array of all peptide sequence encodings. i.e. an array of
#' peptide sequence encodings.
#'
#' @examples
#'
#' examples <- rhoData$data
#'
#' encoded_examples <- encode_onehot(examples)
#'
#' @export
#' @import assertthat
encode_onehot <- function(AA_seqs){

  assert_that(is.list(AA_seqs))
  assert_that(all(unlist(lapply(AA_seqs, nchar))==nchar(AA_seqs[[1]])),
              msg="Aligned sequeces must have the same length.")


  # Create list of amino acids (ordered alphabetically) to index the one-hot
  # representation.
  num_AAs <- 20
  Amino_acids <- unlist(strsplit("ACDEFGHIKLMNPQRSTVWY", ''))
  AAs_to_num <- as.list(1:20)
  names(AAs_to_num) = Amino_acids

  # Convert to character matrix representation.
  AA_matrix <- strs_to_matrix(AA_seqs)

  num_seqs <- length(AA_seqs)
  num_residues <- nchar(AA_seqs[[1]])

  one_hot_array <- array(0L, dim= c(num_seqs, num_residues, num_AAs))

  # For each amino acid in each peptide, add a one at the appropriate position
  # in the one-hot representation.
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
#' @param AA_seqs  A list of aligned peptides.Represented as a list of character
#' vectors.
#'
#' @return Rank 3 array containing all peptide sequence encodings (rank 2).
#' i.e. an array of peptide sequence encodings.
#'
#' @examples
#'
#' examples <- rhoData$data
#'
#' encoded_examples <- encode_physchem(examples)
#' encoded_examples[1]
#'
#' @export
#' @import assertthat
encode_physchem <- function(AA_seqs){
  num_vhses <- 8

  get_VHSEs <- function(AA){
    res <- numeric(num_vhses)
    for (v_num in seq_along(res)){
      res[v_num] <- VHSE[[v_num]][[AA]]
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
