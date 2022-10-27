


#' @importFrom stats model.matrix na.fail prcomp
.run_mat_pca <- function(X, partial_pca=FALSE, initial_dims=50,
                         pca_center=TRUE, pca_scale=FALSE, verbose = NULL) {

  if(verbose) cat("Performing PCA\n")
  if(partial_pca){
    if (!requireNamespace("irlba", quietly = TRUE)) {stop("Package \"irlba\" is required for partial PCA. Please install it.", call. = FALSE)}
    X <- irlba::prcomp_irlba(X, n = initial_dims, center = pca_center, scale = pca_scale)$x

  }else{
    if(verbose & min(dim(X))>2500) cat("Consider setting partial_pca=TRUE for large matrices\n")
    X <- prcomp(X, retx=TRUE, center = pca_center, scale. = pca_scale, rank. = initial_dims)$x
  }

  return(x)

}


.get_shared_mnns <- function(X_mat, Y_mat, XY_shared, nk1 = 10, nk2 = 10) {

  mnns <- batchelor::findMutualNN(XY_shared[[2]], XY_shared[[1]], nk1, nk2)
  mnn_mat <- as.data.frame(mnns)

  Z1 <- sapply(seq_len(nrow(Y_mat)),function(i) mnn_mat[(mnn_mat[,1]==i),2])  #nns of data2 in data1
  Z2 <- sapply(seq_len(nrow(X_mat)),function(i) mnn_mat[(mnn_mat[,2]==i),1])  #nns of data1 in data2

  YX_mnns<-do.call(rbind,lapply(Z1,
                                function(x)
                                  if(length(x)==0) rep(0, nk1)
                                else if(length(x) < nk1) c(x,rep(0, nk1-length(x)))
                                else x)) # nns of data 2 in data1

  XY_mnns<-do.call(rbind,lapply(Z2,
                                function(x)
                                  if(length(x)==0) rep(0, nk2)
                                else if(length(x) < nk2) c(x, rep(0, nk2-length(x)))
                                else x))
  XY_mnns[XY_mnns==0]<-NaN
  YX_mnns[YX_mnns==0]<-NaN

  return(list(XY = XY_mnns, YX = YX_mnns))

}


.check_args <- function(X, Y, coupled_period=30, uncoupled_period=20,
                        dims = 2, perplexity=30, theta=0.5, max_iter = 1000,
                        partial_pca=FALSE, initial_dims=50,
                        pca_center=TRUE, pca_scale=FALSE, Y_init=NULL,
                        stop_lying_iter=ifelse(is.null(Y_init),250L,0L),
                        mom_switch_iter=ifelse(is.null(Y_init),250L,0L),
                        momentum=0.5, final_momentum=0.8, binding_force=0.05,
                        eta=200.0, exaggeration_factor=12.0, num_threads=1,
                        verbose = TRUE) {

  if (!is.matrix(X)) { stop("Input X is not a matrix")}
  if (!is.matrix(Y)) { stop("Input X is not a matrix")}
  if (!(is.logical(pca_center) && is.logical(pca_scale)) ) { stop("pca_center and pca_scale should be TRUE or FALSE")}
  if (!is.wholenumber(initial_dims) || initial_dims<=0) { stop("Incorrect initial dimensionality.")}

  tsne.args <- .check_tsne_params(nrow(X), dims=dims, perplexity=perplexity, theta=theta, max_iter=max_iter, verbose=verbose,
                                  Y_init=Y_init, stop_lying_iter=stop_lying_iter, mom_switch_iter=mom_switch_iter,
                                  momentum=momentum, final_momentum=final_momentum, eta=eta, exaggeration_factor=exaggeration_factor,
                                  binding_force=binding_force,coupled_period=coupled_period,uncoupled_period=uncoupled_period)

  return(tsne.args)

}


#' @export
prepareSATSNE <- function(x_sce, y_sce, shared_feats, assay_type = "logcounts",
                          dimred = NULL, pca = is.null(dimred), shared_pca = TRUE,
                          cosine_norm = TRUE, n_dimred = 50,
                          xy_mnns = NULL, yx_mnns = NULL, nk1 = 10, nk2 = 10,
                          coupled_period=30, uncoupled_period=20, normalise = TRUE,
                          num_threads =1, verbose = TRUE, return_logs = FALSE, ...) {

  # Extract individual matrices
  X_mat <- .get_mat_from_sce(x_sce, exprs_values=assay_type, dimred=dimred, n_dimred=n_dimred)
  Y_mat <- .get_mat_from_sce(y_sce, exprs_values=assay_type, dimred=dimred, n_dimred=n_dimred)


  # Check TSNE arguments are valid
  if(verbose) print("Checking TSNE arguments...")
  tsne_args <- .check_args(X_mat, Y_mat,
                           coupled_period, uncoupled_period, ...)


  # Perform PCA on individual matrices e.g. if dimred = NULL
  if(pca) {
    if(verbose) print("Running PCA on individual matrices.")
    X_mat <- .run_mat_pca(X_mat, ...)
    Y_mat <- .run_mat_pca(Y_mat, ...)
  }


  # TODO: Add checks to make sure shared_feats 2D and correspond to X,Y features
  # Extract shared matrices
  if(verbose) print("Extracting shared feature matrices...")
  X_shared <- .get_mat_from_sce(x_sce[as.character(shared_feats[,1]),],
                                exprs_values=assay_type, dimred=NULL, n_dimred=n_dimred)

  Y_shared <- .get_mat_from_sce(y_sce[as.character(shared_feats[,2]),],
                                exprs_values=assay_type, dimred=NULL, n_dimred=n_dimred)


  # Make feature names consistent
  rownames(Y_shared) <- rownames(X_shared)


  # Apply cosine normalisation
  if(cosine_norm) {
    if(verbose) print("Performing cosine normalisation on the shared matrices.")

    X_l2 <- cosineNorm(X_shared, mode="l2norm")
    X_shared <- .apply_cosine_norm(X_shared, X_l2)

    Y_l2 <- cosineNorm(Y_shared, mode="l2norm")
    Y_shared <- .apply_cosine_norm(Y_shared, Y_l2)
  }


  # Perform PCA jointly on shared feature space
  if(verbose) print("Performing PCA jointly across both datasets on the shared feature space.")
  XY_shared <- batchelor::multiBatchPCA(X_shared, Y_shared, assay.type = assay_type)


  # Compute MNNs in shared feature space
  if(is.null(xy_mnns) && is.null(yx_mnns)) {
    if(verbose) print("Computing MNNs from the shared feature space.")
    mnns <- .get_shared_mnns(X_mat, Y_mat, XY_shared, nk1 = nk1, nk2 = nk2)
    XY_mnn <- mnns$XY
    YX_mnn <- mnns$YX
  }


  # Zero-center individual matrices
  if (normalise) {
    if(verbose) print("Zero-centering the individual matrices.")
    X_mat <- normalize_input(X_mat)
    Y_mat <- normalize_input(Y_mat)
  }


  # Transpose matrices for rapid column-major access.
  if(verbose) print("Transposing matrices for rapid column access.")
  X_mat <- t(X_mat)
  Y_mat <- t(Y_mat)

  return(list(X_mat = X_mat, Y_mat = Y_mat,
              XY_mnn = XY_mnn, YX_mnn = YX_mnn,
              tsne_args = tsne_args))

}


#' TODO: Add BPPARAMS
#' @import Rcpp
#' @export
runSATSNE <- function(x_sce, y_sce, shared_feats, assay_type = "logcounts",
                      dimred = NULL, pca = is.null(dimred), shared_pca = TRUE,
                      cosine_norm = TRUE, n_dimred = 50,
                      xy_mnns = NULL, yx_mnns = NULL, nk1 = 10, nk2 = 10,
                      coupled_period=30, uncoupled_period=20, normalise = TRUE,
                      num_threads =1, verbose = TRUE, return_logs = FALSE, ...) {


  args <- prepareSATSNE(x_sce, y_sce, shared_feats, assay_type = assay_type,
                        dimred = dimred, pca = pca, shared_pca = shared_pca,
                        cosine_norm = cosine_norm, n_dimred = n_dimred,
                        xy_mnns = xy_mnns, yx_mnns = yx_mnns, nk1 = nk1, nk2 = nk2,
                        coupled_period=coupled_period, uncoupled_period=uncoupled_period,
                        normalise = normalise, num_threads =num_threads,
                        verbose = verbose, return_logs = return_logs, ...)

  # Run SATSNE
  if(verbose) print("Starting SATSNE algorithm...")
  out <- do.call(Rsatsne_cpp,
                 c(list(X1=args$X_mat, X2=args$Y_mat, args$XY_mnn, args$YX_mnn,
                        nk1, nk2, distance_precomputed = FALSE,
                        num_threads=num_threads, return_logs = return_logs),
                   args$tsne_args))

  names(out) <- c("Y1", "Y2", "costs1", "costs2", "itercosts", "logs")

  # Transposing back to cells x dims.
  out$Y1 <- t(out$Y1)
  out$Y2 <- t(out$Y2)
  names(out) <- NULL

  return(out)

}


