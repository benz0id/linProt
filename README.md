
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
library(linProt)
#> 
#> Attaching package: 'linProt'
#> The following object is masked from 'package:stats':
#> 
#>     predict
data(rhoData)

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


# View the expected influence of each residue on the function, (lambda max
# in this case).

residue_effect_heatmap(model, 380, 410)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot_cost_over_rep(model)
#> Warning: Ignoring unknown aesthetics: line
#> Ignoring unknown aesthetics: line

# Make some predictions.
results <- predict(valid_data, model)
results[1:5]
#> [1] 546.8277 537.5129 534.1281 541.4704 545.6424
valid_labels[1:5]
#> [1] 564 495 523 551 553
```
