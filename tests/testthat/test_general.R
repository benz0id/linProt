set.seed(1006322089)

test_that("onehot logic functional", {
  seq <- rhoData$data[1]
  one_hot <- encode_onehot(seq)
  num_gaps <- sum(strs_to_matrix(seq) == '-')
  num_residues <- nchar(seq[[1]]) - num_gaps

  expect_equal(sum(one_hot), num_residues)
})

test_that("VHSE correct dimensionality", {
  seq <- rhoData$data[1]
  encoding <- encode_physchem(seq)

  expect_equal(dim(encoding), c(1, nchar(seq[[1]]), 8))
})

test_that("VHSE occupied", {
  seq <- rhoData$data[1]
  encoding <- encode_physchem(seq)
  num_gaps <- sum(strs_to_matrix(seq) == '-')

  expect_true(sum(encoding != 0) > 0.5 * (nchar(seq[[1]]) - num_gaps) * 8)
})


test_that("model does not diverge", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                          encode=encode_onehot)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                       train_labels,
                       valid_data,
                       valid_labels)

  expect_lt(2 * min(model$valid_loss), max(model$valid_loss))
})

test_that("ridge runs", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_physchem)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                        train_labels,
                        valid_data,
                        valid_labels,
                        reg = 'ridge')
  expect_lt(2 * min(model$valid_loss), max(model$valid_loss))
})

test_that("lasso runs", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_physchem)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                        train_labels,
                        valid_data,
                        valid_labels,
                        reg = 'lasso',
                        learning_rate = 0.001)
  expect_lt(2 * min(model$valid_loss), max(model$valid_loss))
})

test_that("elastic runs", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_physchem)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                        train_labels,
                        valid_data,
                        valid_labels,
                        reg = 'lasso',
                        learning_rate = 0.001)
  expect_lt(2 * min(model$valid_loss), max(model$valid_loss))
})

test_that("model can diverge, gradients correct", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_onehot)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  expect_error(linear_train(train_data,
                            train_labels,
                            valid_data,
                            valid_labels,
                            learning_rate = 1))

})

test_that("predictions 5 times better than random for quick model", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_physchem)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                        train_labels,
                        valid_data,
                        valid_labels)
  best_loss <- min(model$valid_losses)
  max_lab <- max(valid_labels)
  min_lab <- min(valid_labels)
  dif <- max_lab - min_lab
  rand <- min_lab + dif * runif(length(valid_labels))

  rand_mean_sq <- mse(rand, valid_labels)

  expect_gt(rand_mean_sq, 5 * best_loss)

})

test_that("predictions 10 times better than random for better model", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_physchem)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                        train_labels,
                        valid_data,
                        valid_labels,
                        learning_rate = 0.01,
                        num_iter = 3000)
  best_loss <- min(model$valid_losses)
  max_lab <- max(valid_labels)
  min_lab <- min(valid_labels)
  dif <- max_lab - min_lab
  rand <- min_lab + dif * runif(length(valid_labels))

  rand_mean_sq <- mse(rand, valid_labels)

  expect_gt(rand_mean_sq, 10 * best_loss)

})

test_that("new sequence is maximal", {
  examples <- rhoData$data
  labels <- rhoData$labels

  # Shuffle, partition, and encode data set.
  shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                           encode=encode_onehot)
  train_data <- shuffled_datasets$e1
  train_labels <- shuffled_datasets$l1
  valid_data <- shuffled_datasets$e2
  valid_labels <- shuffled_datasets$l2

  # Train a linear model to perform regression.
  model <- linear_train(train_data,
                        train_labels,
                        valid_data,
                        valid_labels,
                        learning_rate = 0.02,
                        num_iter = 1000)

  # This is not a great sequence, but it is maximal.
  max_seq <- maximal_sequence(model)
  encoded_max_seq <- encode_onehot(list(max_seq))

  func <- predictions(encoded_max_seq, model)

  expect_gt(func, max(rhoData$labels))

  # Again, this is not a great sequence, but it is minimal
  min_seq <- maximal_sequence(model, do_min = TRUE)
  encoded_min_seq <- encode_onehot(list(min_seq))

  func <- predictions(encoded_min_seq, model)

  expect_lt(func,  min(rhoData$labels))
})

# [END]



