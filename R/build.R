#' Produce a peptide with maximimal predicted function.
#'
#' Uses the weights leaned by a model to produce an amino acid sequence that
#' produces the maximal predicted protein function. Will only work on alignments
#' with low-indel content.
#'
#' Note that a random base is chose when there is a tie between weights.
#'
#' @param model The model as learned by linear_train. This model must have
#' been trained on one-hot encoding.
#'
#' @param do_min Whether to change behaviour to instead minimize the protein
#' function.
#'
#' @return The amino acid sequence of the peptide that maximizes function.
#'
#' @export
maximal_sequence <- function(model, do_min = FALSE){

  w <- model$weights
  assert(dim(w)[2] == 20)

  # Choose to either perform maximization or minimization.
  if (do_min){
    fxn <- min
  } else{
    fxn <- max
  }

  Amino_acids <- unlist(strsplit("ACDEFGHIKLMNPQRSTVWY", ''))

  # Instantiate empty vector to contain final sequence..
  AA_seq <- character(dim(w)[1])

  # Iterate through the weights, choosing the amino acid the either maximizes
  # or minimizes the predicted function of the sequence.
  for (res in seq(dim(w)[1])){
    aa_weights <- w[res, ]

    max_aa <- which(aa_weights == fxn(aa_weights))

    if(length(max_aa) > 1){
      max_aa <- sample(max_aa, 1)
    }
    AA_seq[res] <- Amino_acids[max_aa]
  }

  return(paste(AA_seq, collapse=''))
}
