
<!-- README.md is generated from README.Rmd. Please edit that file -->

# linProt

<!-- badges: start -->
<!-- badges: end -->

linProtFunc provides a pipeline that allows for that seamless training
and evaluation of linear models for protein function prediction from
user provided labelled AA sequences, including support for regression
tasks. This package provides visualizations of the patterns learned by
the linear models which can be used to generate novel proteins with
optimised functional characteristics.

## Installation

You can install the development version of linProt like so:

``` r
install_github("benz0id/linProt")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
data(rhoData)

examples <- rhoData$data
labels <- rhoData$labels

# Shuffle, partition, and encode data set.
shuffled_datasets <- shuffled_partitions(examples, labels, 650,
                                        encode=encode_onehot)

# Train a linear model to perform regression.
model <- linear_train(shuffled_datasets$e1,
                      shuffled_datasets$l1,
                      shuffled_datasets$e2,
                      shuffled_datasets$l2)

# View the expected influence of each residue on the function, (lambda max
# in this case).

residue_effect_heatmap(model, 380, 410)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot_cost_over_rep(model)

# Make some predictions.
results <- predict(shuffled_datasets$e2, model)
results[1:5]
#> [1] 540.1165 525.5414 529.7755 538.2172 532.9388
shuffled_datasets$l2[1:5]
#> [1] 536 540 530 544 533
```
