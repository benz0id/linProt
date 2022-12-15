
<!-- README.md is generated from README.Rmd. Please edit that file -->

# linProt

<!-- badges: start -->
<!-- badges: end -->

## Description

This package provides a pipeline that allows for that training and
evaluation of linear models for protein function prediction from user
provided labelled AA sequences. This package provides visualizations of
the patterns learned by the linear models and provides functions that
can be used to generate novel proteins with optimised functional
characteristics.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("benz0id/linProt", build_vignettes = TRUE)
library("linProt")
```

## Overview

``` r
ls("package:linProt")
data(package = "linProt")
browseVignettes("linProt")
```

This package works with continuously distributed protein functional
characteristics, such as peak light absorption wavelength or
fluorescence.

In the context of protein function prediction, linear models assign
weights to each of the amino acid positions along an alignment. We can
interpret these weights to gain insight into how residues contribute to
the protein’s function and use the model to make predictions on unseen
sequences.

This model takes a list of aligned amino acid sequences, each with a
label denoting their experimentally determined function. These sequences
are then encoded as either as A x 20 one hot matricies (where A is the
number of residues in the aligned protein) or as A x 8 matricies, using
the VHSE8 scores of each residue. After being partitioned into training
and validation sets, these encoded data are then used to train a
regularised linear model.

This trained linear model can be used for the function prediction of
novel sequences, the design of maximally or minimally functional
proteins, or to gain insight into the functional characteristics of the
protein.

In order to achieve this, the package implements 8 functions that can be
ordered form a modular pipeline. These include:

***shuffled_partitions*** - Shuffle, partition, and encode the supplied
peptide alignment using one of the following encoding functions

***encode_onehot*** - Encode a peptide alignment as a tensor of one-hot
vectors.

***encode_physchem*** - Encode a peptide alignment as a tensor of
vectors encoded by the VHSE matrix.

***linear_train*** - Train and return a linear model that predicts
protein function on using some training and validation data given some
hyperparameters.

***predictions*** - Use a linear model to make function predictions for
some aligned peptide sequences.

***maximal_sequence*** - Create a peptide that would have the maximum
(or minimum) possible function as predicted by some model.

***plot_cost_over_rep*** - Show the change in cost over the duration of
some linear model’s training cycles. Helpful for discovering
overfitting.

***property_effect_heatmap*** - Using the weights supplied by a model
trained on VHSE8 encoded peptides, show the learned effect of each VHSE
on protein function.

***residue_effect_heatmap*** - Using the weights supplied by a model
trained on one hot encoded peptides, show the learned effects that
having any given residue at any given position will have on protein
function.

![](./inst/extdata/diagram.png)

## Example

With this package, you can easily generate, evaluate, and use linear
models for protein function/fitness predictionn tasks.

Here, we generate a linear model to peodict the absorption properties of
various aligned channel rhodopin proteins using the provided dataset.

``` r
library(linProt)
data(rhoData)

examples <- rhoData$data
labels <- rhoData$labels

examples[1]
#> $`1`
#> [1] "--------------------------------------------------mfai----npeymn---------------------------et----vlld-------ec--tpiy-----ldigplweqvvar--vtqwfg-vilslvfliyyiwntykatcg-weelyvctvefckiiielyf-eytppam-------------------ifqtngqvtpwlryaewlltcpvilihls-nitgl--------ndd---ysgrtms-litsdlggicmavtaalsk-------g-----wlka-lffvigcg-ygastfynaaciyiesy----------------y------------------tm--p---qgicrrlvlwmagvfftswfmfpglflagpe-gtq---------al-swagttightvadllsknawgmighflrveihkhiiihgd-------vrrpvtvkalgrqvsvncfvdke------eeeederi--------"
labels[1]
#> [1] 436

# Shuffle, partition, and encode data set.
shuffled_datasets <- shuffled_partitions(examples, labels, num_p1 = 650,
                                        encode=encode_onehot)

# Train a linear model to perform regression.
model <- linear_train(train_data = shuffled_datasets$e1,
                      train_labels = shuffled_datasets$l1,
                      valid_data = shuffled_datasets$e2,
                      valid_labels = shuffled_datasets$l2)

# View the expected influence of each residue on the function, (lambda max
# in this case).

generated_plot <- residue_effect_heatmap(model, 380, 410)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
generated_plot <- plot_cost_over_rep(model)

# Make some predictions.
results <- predictions(shuffled_datasets$e2, model)
as.integer(results[1:5])
#> [1] 541 537 543 469 414
shuffled_datasets$l2[1:5]
#> [1] 526 495 573 470 515
```

See the vignettes for a detailed description of the available functions,
which can provide a manner of different visualizations of the model’s
attributes, as well as apply it to protein engineering tasks.

## Contributions

This package was written by Benjamin Tudor Price. The R package
`ggplot2` was used for visualization of gradient descent progress in the
*plot_cost_over_rep* function. The R package `gplots` was used to
visualise the learned weights using a heatmap for the
*property_effect_heatmap* and *residue_effect_heatmap* functions. The R
packages `assertthat` and `assert` were used to check user input
throughout the project. In the function *linear_train*, The R package
`progess` was used to produce an informative progress bar. The R package
`stats` was also used in *linear_train* in order to help format
variables. The `shiny` package was used to develop the shiny app.

The included dataset used algined sequence data published by (Karasuyama
et al., 2018)

The VHSE matrix used (VHSE8) was originally released by (Mei et al.,
2005)

## References

Mei, H., Liao, Z. H., Zhou, Y., & Li, S. Z. (2005). A new set of amino
acid descriptors and its application in peptide QSARs. Biopolymers,
80(6), 775–786. <https://doi.org/10.1002/bip.20296>

Xu, Y., Verma, D., Sheridan, R. P., Liaw, A., Ma, J., Marshall, N. M.,
McIntosh, J., Sherer, E. C., Svetnik, V., & Johnston, J. M. (2020). Deep
Dive into Machine Learning Models for Protein Engineering. Journal of
Chemical Information and Modeling, 60(6), 2773–2790.
<https://doi.org/10.1021/acs.jcim.0c00073>

Karasuyama, M., Inoue, K., Nakamura, R., Kandori, H., & Takeuchi, I.
(2018). Understanding Colour Tuning Rules and Predicting Absorption
Wavelengths of Microbial Rhodopsins by Data-Driven Machine-Learning
Approach. Scientific Reports, 8(1).
<https://doi.org/10.1038/s41598-018-33984-w>

Bedbrook, C. N., Yang, K. K., Robinson, J. E., Mackey, E. D., Gradinaru,
V., & Arnold, F. H. (2019). Machine learning-guided channelrhodopsin
engineering enables minimally invasive optogenetics. Nature Methods,
16(11), 1176–1184. ://doi.org/10.1038/s41592-019-0583-8

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
New York, 2016.

Warnes G, Bolker B, Bonebakker L, Gentleman R, Huber W, Liaw A, Lumley
T, Maechler M, Magnusson A, Moeller S, Schwartz M, Venables B (2022).
*gplots: Various R Programming Tools for Plotting Data*. R package
version 3.1.3, <https://CRAN.R-project.org/package=gplots>.

Wickham H (2019). *assertthat: Easy Pre and Post Assertions*. R package
version 0.2.1, <https://CRAN.R-project.org/package=assertthat>.

Binette O (2020). *assert: Validate Function Arguments*. R package
version 1.0.1, <https://CRAN.R-project.org/package=assert>.

Csárdi G, FitzJohn R (2019). *progress: Terminal Progress Bars*. R
package version 1.2.2, <https://CRAN.R-project.org/package=progress>.

R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria. URL
<https://www.R-project.org/>.

Yihui Xie (2022). knitr: A General-Purpose Package for Dynamic Report
Generation in R. R package version 1.40.

## Ackowledgements

This package was developed as part of an assessment for 2022 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `linProt` welcomes issues, enhancement requests, and other
contributions. This package was improved by the advice of my peers and
Dr. Anjali Silva. To submit an issue, use the GitHub issues.
