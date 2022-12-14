#' Heatmap View of Learned One Hot Weights
#'
#'
#' A function that generates an N x 20 heatmap displaying the weight assigned to
#' each possible residue across a specified region of the alignment. Allows for
#' the  interpretation of the learned attributes of the dataset. e.g. a red tile
#' at a given position indicates that amino acid is expected to increase the
#' protein function when present at that position.
#'
#' Extend the width of the produced figure as necessary for best results.
#'
#' @param model The model as learned by linear_train. This model must have
#' been trained on a one-hot oncoded training set.
#'
#' @param start The first in the range of residues to be displayed, 1 by
#' default. Must be less than start.
#'
#' @param stop The last in the range of residues to be displayed, the last
#' residue in the alignment by default. Must be greater than start. If left
#' unspecified, defaults to the full length of the alignment.
#'
#' @examples
#'
#'
#' examples <- rhoData$data
#' labels <- rhoData$labels
#'
#' shuffled_datasets <- shuffled_partitions(examples, labels, 650,
#'                                          encode=encode_onehot)
#'
#' # Train a linear model to perform regression.
#'
#'
#' model <- linear_train(shuffled_datasets$e1,
#'                       shuffled_datasets$l1,
#'                       shuffled_datasets$e2,
#'                       shuffled_datasets$l2)
#'
#'
#' residue_effect_heatmap(model)
#'
#' @export
#' @import gplots
residue_effect_heatmap <- function(model,
                                   start=1,
                                   stop=-1){
  w <- t(model$weights)

  # If stop was unspecified, set it to the maximum length.
  if (stop == -1){
    stop <- dim(w)[2]
  }

  # Check that weights are one-hot encoded.
  assert(dim(w)[1] == 20)

  # Prepare data for visualisation.
  Amino_acids <- unlist(strsplit("ACDEFGHIKLMNPQRSTVWY", ''))
  sel_w <- w[, start:stop]
  colnames(sel_w) <- start:stop
  rownames(sel_w) <- Amino_acids

  heatmap.2(sel_w, scale = "none", col = bluered(100),
            trace = "none", density.info = "none", Colv = FALSE, Rowv = FALSE,
            dendrogram="none",
            xlab="Aligned Residue #", ylab="Amino Acid", main="Effect of Amino
            Acid prescence on Function", key.xlab = '', key.ylab = '')

  return(invisible(NULL))

}


#' Heatmap View of Physiochemical Weighings
#'
#'
#' A function that generates an N x 8 heat map displaying the weight assigned to
#' each of the VHSE8 properties across a specified region of the alignment.
#' Allows for the interpretation of the learned attributes of the dataset. e.g.
#' a red tile at a given position indicates that physiochemical property is
#' related to an increase the protein function when present at that position.
#'
#' Extend produced figure as necessary for best results.
#'
#' @param model The model as learned by linear_train. This model must have
#' been trained on a physiochemical encoding.
#'
#' @param start The first in the range of residues to be displayed, 1 by
#' default. Must be less than start.
#'
#' @param stop The last in the range of residues to be displayed, the last
#' residue in the alignment by default. Must be greater than start. If left
#' unspecified, defaults to the full length of the alignment.
#'
#' @examples
#'
#' examples <- rhoData$data
#' labels <- rhoData$labels
#'
#' shuffled_datasets <- shuffled_partitions(examples, labels, 650,
#'                                          encode=encode_physchem)
#'
#' # Train a linear model to perform regression.
#' model <- linear_train(shuffled_datasets$e1,
#'                       shuffled_datasets$l1,
#'                       shuffled_datasets$e2,
#'                       shuffled_datasets$l2)
#'
#'
#' property_effect_heatmap(model)
#'
#' @export
#' @import gplots
property_effect_heatmap <- function(model,
                                    start=1,
                                    stop=-1){

  w <- t(model$weights)
  assert(dim(w)[1] == 8)

  if (stop == -1){
    stop <- dim(w)[2]
  }

  vhses <- c(
    "VHSE1",
    "VHSE2",
    "VHSE3",
    "VHSE4",
    "VHSE5",
    "VHSE6",
    "VHSE7",
    "VHSE8"
  )
  sel_w <- w[, start:stop]
  colnames(sel_w) <- start:stop
  rownames(sel_w) <- vhses

  heatmap.2(sel_w, col = bluered(100), scale='none',
          trace = "none", density.info = "none", Colv = FALSE, Rowv = FALSE,
          dendrogram='none', key.title = "Relative Weight", keysize=1.5,
          xlab="Residue", ylab="Amino Acid", main="Relative Effect of Each
          VHSE on Function",
          offsetRow = 0,
          offsetCol = 0, key.xlab = '', key.ylab = '')

  return(invisible(NULL))
}


#' Plot the Change in Cost Over Gradient Descent Iterations
#'
#' Visualizes the change in cost over both the training and validation sets
#' throughout the course of the training. Allows for tuning of hyperparameters
#' and can help prevent overfitting.
#'
#' @param model A model produced by linear train.
#'
#' @examples
#' examples <- rhoData$data
#' labels <- rhoData$labels
#'
#' shuffled_datasets <- shuffled_partitions(examples, labels, 650,
#'                                          encode=encode_physchem)
#'
#' # Train a linear model to perform regression.
#'
#'
#' model <- linear_train(shuffled_datasets$e1,
#'                       shuffled_datasets$l1,
#'                       shuffled_datasets$e2,
#'                       shuffled_datasets$l2)
#'
#'
#' plot_cost_over_rep(model)
#'
#' @export
#' @import ggplot2
plot_cost_over_rep <- function(model){

  x <- seq_along(model$train_losses) * model$every

  progress_data <- data.frame(iteration=x,
                              train_cost = log(model$train_losses),
                              valid_cost = log(model$valid_losses))

  ggplot(progress_data, aes(x=iteration, color=dataset)) +
    geom_line(aes(x=iteration, y=valid_cost, colour='validation')) +
    geom_line(aes(x=iteration, y=train_cost, colour='training')) +
    labs(title="Change in Cost Through Gradient Descent",
         x="Iteration", y = "Log Cost")
}

# [END]
