
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
# install.packages("devtools")
devtools::install_github("benz0id/linProt", build_vignettes = TRUE)
library("benz0id/linProt")
```

## Example

With this package, you can easily generate, evaluate, and use linear
models for protein function/fitness predictionn tasks.

Here, we generate a linear model to peodict the absorption properties of
various aligned channel rhodopin proteins using the provided dataset.

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

examples[1:5]
#> $`1`
#> [1] "--------------------------------------------------mfai----npeymn---------------------------et----vlld-------ec--tpiy-----ldigplweqvvar--vtqwfg-vilslvfliyyiwntykatcg-weelyvctvefckiiielyf-eytppam-------------------ifqtngqvtpwlryaewlltcpvilihls-nitgl--------ndd---ysgrtms-litsdlggicmavtaalsk-------g-----wlka-lffvigcg-ygastfynaaciyiesy----------------y------------------tm--p---qgicrrlvlwmagvfftswfmfpglflagpe-gtq---------al-swagttightvadllsknawgmighflrveihkhiiihgd-------vrrpvtvkalgrqvsvncfvdke------eeeederi--------"
#> 
#> $`2`
#> [1] "----------------------------------------------------------------------------------------------------------------------------gtngaqtasn--vlqwqlaagfsilllmfyayqtwkstcg-weeiyvcaiemvkvilefff-efknpsm-------------------lylatghrvqwlryaewlltcpvilihls-nltgl--------snd---dssrtmg-llacsigtivwgatsamas-------g-----yvkv-iffclgvy-c-antf-yraqayikgy----------------h------------------tv--p---kgrcrqvvtgmawlffvswgmfpilfilgpe-gfg---------vl-svygstvghtiidpmsknrcglpghyprvlv-------------------------------------------------------------"
#> 
#> $`3`
#> [1] "---------------------------------------------------------------------------------------------------------------------------ngtnaeklaan--ilqwit-falsalclmfygyqtwkstcg-weeiyvatiemikfiieyfh-efdepav-------------------iyssngnktvwlryaewlltcpvilihls-nltgl--------and---ynkrtmg-llvsdigtivwgttaalsk-------g-----yvrv-ifflmglc-ygiytffnaakvyieay----------------h------------------tv--p---kgicrdlvrylawlyfcswamfpvlfllgpe-gfg---------hi-nqfnsaiahaildlasknawsmmghflrv---------------------------------------------------------------"
#> 
#> $`4`
#> [1] "---------------------------------------------------------------------------------------------------------------mdpia---lq--agydllgdgrpetlwlgigtllmligtfyf-lvrgwgvtdkdareyyavtilvpgiasaaylsmffgigltev---------------tvg-gemldiyyaryadwlfttplllldla-llakv--------------drvtigt-lvgvdalaivtvligalsh------ta-----iary-swwlfsti-cmivvlyflatslrsaa----------------k----------------------er--gpevastfntltalvlvlwtaypilwiigte-gag---------vv-glgietllfmvldvttkvgfgfillrsrailgdteapeps------------agadvs------aad-------------------------"
#> 
#> $`5`
#> [1] "-----msrrpwllalalavalaagsagast-----------gsdatvpvatq-d----gpdy----vfhrahermlfqts--ytlen---ngsvicipnng-------qcfclawl-----ksngtnaeklaan--ilqwit-falsalclmfygyqtwkstcg-weeiyvatiemikfiieyfh-efdepav-------------------iyssngnktvwlryaewlltcpvilihls-nltgl--------and---ynkrtmg-llvsdiggivwattaalsk-------g-----yvrv-ifflmglc-ygiytffnaakvyieay----------------h------------------tv--p---kgrcrqvvtgmawlffvswgmfpilfilgpe-gfg---------vl-svygstvghtiidlmskncwgllghylrvlihehilihgd-------irkttklniggteievetlvede------aeagavssedlyfq--"
labels[1:5]
#> [1] 436 444 453 455 455

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

residue_effect_heatmap(model, 380, 410)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot_cost_over_rep(model)
#> Warning: Ignoring unknown aesthetics: line
#> Ignoring unknown aesthetics: line

# Make some predictions.
results <- predict(shuffled_datasets$e2, model)
as.integer(results[1:5])
#> [1] 591 430 536 539 478
shuffled_datasets$l2[1:5]
#> [1] 606 453 516 540 540
```

See the vignettes for a detailed description of the available functions,
which can provide a manner of different visualizations of the model’s
attributes, as well as apply it to protein engineering tasks.
