#' Linprot
#'
#' linProt provides a pipeline that allows for that seamless training and
#' evaluation of linear models for protein function prediction from user
#' provided labelled AA sequences, including support for regression  tasks.
#' This package provides visualizations of the patterns learned by the linear
#' models which can be used to generate novel proteins with optimized
#' functional characteristics.
#'
#' @section linProt functions:
#' The linProt package provides 10 functions:
#' \itemize{
#'   \item shuffled_partitions
#'   \item encode_onehot
#'   \item encode_physchem
#'   \item linear_train
#'   \item predictions
#'   \item maximal_sequence
#'   \item plot_cost_over_rep
#'   \item property_effect_heatmap
#'   \item residue_effect_heatmap
#'   \item runLinProt
#' }
#'
#' @section linProt data:
#' The linProt packages also provides a basic dataset for the training of
#' linear models, "rhoData".
#'
#' @author {Benjamin Tudor Price, \email{bejamin.tudorprice@mail.utoronto.ca}. }
#'
#' @examples
#' data(rhoData)
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
#' # Train a linear model to perform regression.
#' model <- linear_train(train_data,
#'                       train_labels,
#'                       valid_data,
#'                       valid_labels)
#'
#'
#' # View the expected influence of each residue on the function, (lambda max
#' # in this case).
#' residue_effect_heatmap(model, 380, 410)
#'
#' cost_graph <- plot_cost_over_rep(model)
#'
#' results <- predictions(valid_data, model)
#'
#' @references
#' Mei, H., Liao, Z. H., Zhou, Y., & Li, S. Z. (2005). A new set of amino acid
#' descriptors and its application in peptide QSARs. Biopolymers, 80(6),
#' 775–786. https://doi.org/10.1002/bip.20296
#'
#' Xu, Y., Verma, D., Sheridan, R. P., Liaw, A., Ma, J., Marshall, N. M.,
#' McIntosh, J., Sherer, E. C., Svetnik, V., & Johnston, J. M. (2020).
#' Deep Dive into Machine Learning Models for Protein Engineering.
#' Journal of Chemical Information and Modeling, 60(6), 2773–2790.
#' https://doi.org/10.1021/acs.jcim.0c00073
#'
#' Karasuyama, M., Inoue, K., Nakamura, R., Kandori, H., & Takeuchi, I. (2018).
#' Understanding Colour Tuning Rules and Predicting Absorption Wavelengths of
#' Microbial Rhodopsins by Data-Driven Machine-Learning Approach. Scientific
#' Reports, 8(1). https://doi.org/10.1038/s41598-018-33984-w
#'
#' Bedbrook, C. N., Yang, K. K., Robinson, J. E., Mackey, E. D., Gradinaru, V.,
#' & Arnold, F. H. (2019). Machine learning-guided channelrhodopsin engineering
#' enables minimally invasive optogenetics. Nature Methods, 16(11), 1176–1184.
#' https://doi.org/10.1038/s41592-019-0583-8
#'
#' Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R
#' project for statistical computing devoted to biological sequences retrieval
#' and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
#' Structural approaches to sequence evolution: Molecules, networks,
#' opulations, series Biological and Medical Physics, Biomedical Engineering,
#' 207-232. Springer Verlag, New York. ISBN : 978-3-540-35305-8.
#'
#' @docType package
#' @name linProt
NULL
#> NULL
