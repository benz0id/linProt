


#' Produce a Sequence to Optimize a Given Function
#'
#' Uses the weights leaned by a model to produce an amino acid sequence that
#' produces the maximal predited protein function. Will only work on alignments
#' with low-indel content.
#'
#' Note that a random base is chose when there is a tie between weights.
#'
#' @param model The model as learned by linear_train. This model must have
#' been trained on one-hot encoding.
#'
#' @param min Whether to change behaviour to instead minimise the protein
#' function.
#'
#' @return A character representation of the best peptide.
#'
#'
#' @export
maximal_sequence <- function(model, min = FALSE){
  w <- model$weights

  if (min){
    fxn <- min
  } else{
    fxn <- max
  }

  Amino_acids <- unlist(strsplit("ACDEFGHIKLMNPQRSTVWY", ''))

  AA_seq <- character(dim(w)[1])

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
