#' MPLNClust: Clustering via mixtures of multivariate Poisson-log normal distribution
#'
#' linProtFunc provides a pipeline that allows for that seamless training and
#' evaluation of linear models for protein function prediction from user
#' provided labelled AA sequences, including support for regression  tasks.
#' This package provides visualizations of the patterns learned by the linear
#' models which can be used to generate novel proteins with optimised
#' functional characteristics.
#'
#' @section MPLNClust functions:
#' The MPLNClust package provides 10 functions:
#' \itemize{
#'   \item shuffled_partitions
#'   \item encode_onehot
#'   \item encode_physchem
#'   \item linear_train
#'   \item predict
#'   \item maximal_sequence
#'   \item plot_cost_over_rep
#'   \item property_effect_heatmap
#'   \item residue_effect_heatmap
#' }
#'
#' @author {Benjamin Tudor Price, \email{bejamin.tudorprice@mail.utoronto.ca}. }
#'
#' @examples
#'
#' example_data <- get_example_data()
#'
#' examples <- example_data$data
#' labels <- example_data$labels
#'
#' # Shuffle, partition, and encode data set.
#'
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
#'
#' \dontrun{
#' model <- linear_train(train_data,
#'                       train_labels,
#'                       valid_data,
#'                       valid_labels)
#' }
#'
#' # View the expected influence of each residue on the function, (lambda max
#' # in this case).
#'
#' residue_effect_heatmap(model, 380, 410)
#'
#' plot_cost_over_rep(model)
#' @references
#'
#'
#'
#'
#'
#'
#'
#'
#' Silva, A. et al. (2019). A multivariate Poisson-log normal mixture model
#' for clustering transcriptome sequencing data. \emph{BMC Bioinformatics} 20.
#' \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2916-0}{Link}
#'
#' Subedi, S., and R. Browne (2020). A parsimonious family of multivariate Poisson-lognormal
#' distributions for clustering multivariate count data. arXiv preprint arXiv:2004.06857.
#' \href{https://arxiv.org/pdf/2004.06857.pdf}{Link}
#'
#' Aitchison, J. and C. H. Ho (1989). The multivariate Poisson-log normal distribution.
#' \emph{Biometrika} 76.
#'
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#' principle. In \emph{Second International Symposium on Information Theory}, New York, NY,
#' USA, pp. 267–281. Springer Verlag.
#'
#' Arlot, S., Brault, V., Baudry, J., Maugis, C., and Michel, B. (2016).
#' capushe: CAlibrating Penalities Using Slope HEuristics. R package version 1.1.1.
#'
#' Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#' clustering with the integrated classification likelihood. \emph{IEEE Transactions
#' on Pattern Analysis and Machine Intelligence} 22.
#'
#' Bozdogan, H. (1994). Mixture-model cluster analysis using model selection criteria
#' and a new informational measure of complexity. In \emph{Proceedings of the First US/Japan
#' Conference on the Frontiers of Statistical Modeling: An Informational Approach:
#' Volume 2 Multivariate Statistical Modeling}, pp. 69–113. Dordrecht: Springer Netherlands.
#'
#' Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential
#' expression analysis of RNA-seq data. \emph{Genome Biology} 11, R25.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics}
#' 6.
#'
#' @docType package
#' @name linProt
NULL
#> NULL
