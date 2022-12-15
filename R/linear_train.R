
# Helper function for training and prediction.
reduce_input_dim <- function(examples){
  assert_that(is.array(examples))

  # Convert rank 3 tensor to rank 2 tensor by linearizing 2nd and 3rd into one.
  new_shape <- c(dim(examples)[1], dim(examples)[2] * dim(examples)[3])
  new_examples <- array(examples, dim=new_shape)

  # Add dummy biases column for obtaining intercept.
  new_examples <- cbind(new_examples, 1L)
  return(new_examples)
}

# Helper function for training and prediction
predict_linearised <- function(X, w){
  preds <- X %*% w
  return(preds)
}

#' Predict a Continuously Distributed Protein Function
#'
#' Performs regression using the given, trained model in order to predict
#' some continuously distributed protein function.
#'
#' @param X Encoded protein sequences. A rank 3 tensor containing examples with
#' the same dimensionality as those used to train the model.
#'
#' @param model The model used to perform the prediction.
#'
#' @returns Vector of doubles. Predictions for each input example.
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
#'
#' model <- linear_train(shuffled_datasets$e1,
#'                       shuffled_datasets$l1,
#'                       shuffled_datasets$e2,
#'                       shuffled_datasets$l2)
#'
#' predictions(shuffled_datasets$e1[50:100, ,], model)
#'
#' @export
#' @import assertthat
predictions <- function(X, model){
  # Reduce rank of encodings and weights. Add intercept to weights and dummy
  # biases to examples.
  X_lin <- reduce_input_dim(X)
  w <- model$weights
  dim(w) <- dim(w)[1] * dim(w)[2]
  w <- append(w, model$intercept, length(w))

  # Return the prediction.
  return(predict_linearised(X_lin, w))
}

# Helper function for training.
get_deriv <- function(X, w, t){
  y <- predict_linearised(X, w)
  # Calculate derivatives using MSE as the cost function.
  derivs <- t(X) %*% (y - t)
  return(list(derivs, y))
}

# Helper function for training.
regularise <- function(l1, l2, w){
  reg_vector <- numeric(length(w))

  # Calculate l1 weight decay.
  if (l1 != 0){
    for (i in seq_along(w)){
      weight <- w[i]
      if (weight > 0){
        reg_vector[i] <- l1
      } else if (weight < 0){
        reg_vector[i] <- - l1
      } else {
        reg_vector[i] <- 0
      }
    }
  } else {
    # Do not perform l1 regularization.
  }

  # Calculate l2 weight decay.
  reg_vector <- reg_vector + l2 * w
  # Remove bias adjustment.
  reg_vector[length(reg_vector)] <- 0
  return(reg_vector)
}

# Helper function for training.
mse <- function(y, t){
  sq <- sum((y - t) ** 2)
  return(sq)
}

#' Train a Linear Regression Model for Prediction of Protein Function
#'
#' Trains a linear model using gradient descent. Stores the performance on the
#' training and validations sets throughout the process.
#'
#' @param train_data Rank 3 array of encoded protein sequences.
#'
#' @param train_labels Numeric vector of labels for each of the given training
#' examples.
#'
#' @param valid_data Rank 3 encoded array of aligned protein sequences.
#'
#' @param valid_labels Numeric vector of labels for each of the given validation
#' examples.
#'
#' @param reg c('elastic, ridge, lasso'), elastic by default. Specified the
#' regularization technique to be used. Appropriate regularization
#' hyperparameters must be defined.
#'
#' @param reg_hypers Named numeric vector containing the l1 and l2 constants to
#' be used in regularization. Both 0.01 by default. NOTE: Ensure that l1 and l2
#' hyperparameters required for chosen regularization are used.
#'
#' @param num_iter The integer number of iterations to be performed during gradient
#' descent. 1000 by default. Higher values often lead to more performant models,
#' but can also lead to overfitting.
#'
#' @param rec_loss_every Integer. The program will record the cost on the training set
#' every <rec_loss_every> iterations.
#'
#' @param learning_rate The double learning rate for gradient descent.
#'
#' @return A list containing the trained model and its hyperparameters,
#' consisting of:
#' weights - A matrix containing the learned weights.
#' intercept - The bias of the model.
#' train_losses - The total loss of the model on the training set at each step
#' of gradient descent.
#' valid_losses - The total loss of the model on the validation set at each step
#' of gradient descent.
#' every - Equal to num_iter, used to calculate the number of gradient descent
#' iterations.
#'
#'
#' @examples
#' examples <- rhoData$data
#' labels <- rhoData$labels
#'
#' shuffled_datasets <- shuffled_partitions(examples, labels, 650,
#'                                          encode=encode_physchem)
#'
#' model <- linear_train(train_data = shuffled_datasets$e1,
#'                       train_labels = shuffled_datasets$l1,
#'                       valid_data = shuffled_datasets$e2,
#'                       valid_labels = shuffled_datasets$l2,
#'                       reg = 'ridge',
#'                       reg_hypers = setNames(c(0.001, 0.001), c("l1", "l2")),
#'                       num_iter = 500,
#'                       rec_loss_every = 100,
#'                       learning_rate = 0.001
#'                       )
#'
#' predictions(shuffled_datasets$e1[50:60, ,], model)
#'
#' @export
#' @import assert
#' @import assertthat
#' @import progress
#' @importFrom stats setNames
linear_train <- function(train_data,
                         train_labels,
                         valid_data,
                         valid_labels,
                         reg = 'elastic',
                         reg_hypers = setNames(c(0.01, 0.01), c("l1", "l2")),
                         num_iter = 500,
                         rec_loss_every = 10,
                         learning_rate = 0.01
                         ){
  assert_that(is.array(train_data))
  assert_that(is.vector(train_labels))
  assert_that(is.array(valid_data))
  assert_that(is.vector(valid_labels))

  assert_that(num_iter >= 1)
  assert_that(rec_loss_every >= 1)
  assert_that(learning_rate > 0)

  encoding_match <- all(dim(train_data)[2:3] == dim(valid_data)[2:3])
  assert(encoding_match, msg="The training set and validation set do not have
         the same encoding.")

  proper_label_num <- dim(train_data)[1] == length(train_labels) &
    dim(valid_data)[1] == length(valid_labels)
  assert(proper_label_num, msg="One label must be present for each training and
         validation example.")

  l1 <- 0
  l2 <- 0
  # Ensure that the user has the defined l1 and l2 hyperparameters properly,
  # which depend the the regularisation technique specified.
  if(reg =='elastic'){
    right_hypers <- is.null(reg_hypers[["l1"]]) | is.null(reg_hypers[["l2"]])
    assert(!right_hypers, msg="Elastic regularisation requires both l1
                and l2 regularisation hyperparameters.")
    l1 <- reg_hypers[["l1"]]
    l2 <- reg_hypers[["l2"]]
  } else if(reg =='lasso'){
    right_hypers <- is.null(reg_hypers[["l1"]])
    assert(!right_hypers, msg="Lasso regularisation requires both l1
                and l2 regularisation hyperparameters.")
    l1 <- reg_hypers[["l1"]]
  } else if(reg =='ridge'){
    right_hypers <- is.null(reg_hypers[["l2"]])
    assert(!right_hypers, msg="Ridge regularisation requires both l1
                and l2 regularisation hyperparameters.")
    l2 <- reg_hypers[["l2"]]
  }  else {
    stop("unknown regularisation technique. reg can be elastic, ridge or lasso.
         ")
  }

  number_losses_to_rec <- num_iter %/% rec_loss_every

  w_final_dim <- c(dim(train_data)[2], dim(train_data)[3])

  # Linearize examples and add biases column.
  X_train <- reduce_input_dim(train_data)
  X_valid <- reduce_input_dim(valid_data)
  N_train <- dim(X_train)[1]
  N_valid <- dim(X_valid)[1]

  # Initilize variables to store training progress info.
  M <- dim(X_train)[2]

  train_losses <- numeric(number_losses_to_rec)
  valid_losses <- numeric(number_losses_to_rec)

  w <- numeric(M)

  pb <- progress_bar$new(total = num_iter,
                         format = "  train [:bar] :percent eta: :eta :msg :valid_loss")

  # Perform gradient descent.
  valid_loss <- 0
  for(i in seq(num_iter)){
    # Update progress bar.
    pb$tick(tokens = list(msg="Validation Cost:", valid_loss=valid_loss))

    # Compute derivative and update weights.
    res <- get_deriv(X_train, w, train_labels)
    weight_deriv <- res[[1]]
    y <- res[[2]]

    reg_values <- regularise(l1, l2, w)

    adj <- learning_rate * (weight_deriv + reg_values) / (M + 1)
    w <- w - adj

    if (NaN %in% w){
      stop('Weights have diverged. Try more conservative hyperparameters.')
    }
    # Record losses every few iterations.
    if (i %% rec_loss_every == 0){
      yv <- predict_linearised(X_valid, w)
      train_losses[i %/% rec_loss_every] <- (mse(y, train_labels) /
                                               length(train_labels))
      valid_loss <- mse(yv, valid_labels) / length(valid_labels)
      valid_losses[i %/% rec_loss_every] <- valid_loss
    }
  }

  # Store and return model, metadata, and training progress data.
  intercept <- w[length(w)]
  w <- w[-length(w)]
  dim(w) <- w_final_dim
  rtrn <- list(w, intercept, train_losses, valid_losses, rec_loss_every)
  names(rtrn) <- c('weights', 'intercept', 'train_losses', 'valid_losses',
                   'every')
  return(rtrn)
}

# [END]
