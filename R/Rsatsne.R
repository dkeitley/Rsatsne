#' Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
#'
#' Wrapper for the C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding. t-SNE is a method for constructing a low dimensional embedding of high-dimensional data, distances or similarities. Exact t-SNE can be computed by setting theta=0.0.
#'

#' @useDynLib Rtsne, .registration = TRUE
#' @import Rcpp
#' @importFrom stats model.matrix na.fail prcomp
#'
#' @export
Rsatsne <- function (X, ...) {
  UseMethod("Rsatsne", X)
}

#' @describeIn Rtsne Default Interface
#' @export
Rsatsne.default <- function(X,Y, Xshared,Yshared,nk1=10,nk2=10, mat12=NULL,mat21=NULL,dims=2, initial_dims=50,
                            perplexity=30, theta=0.5,
                            check_duplicates=TRUE,
                            pca=TRUE, partial_pca=FALSE, max_iter=1000,verbose=getOption("verbose", FALSE),
                            is_distance=FALSE, Y_init=NULL,
                            pca_center=TRUE, pca_scale=FALSE, normalize=TRUE,
                            stop_lying_iter=ifelse(is.null(Y_init),250L,0L),
                            mom_switch_iter=ifelse(is.null(Y_init),250L,0L),
                            momentum=0.5, final_momentum=0.8, binding_force=0.05,
                            coupled_period=30, uncoupled_period=20,
                            eta=200.0, exaggeration_factor=12.0, num_threads=1, ...) {

  print("Entered Rsatsne...")

  if (!is.logical(is_distance)) { stop("is_distance should be a logical variable")}
  if (!is.matrix(X)) { stop("Input X is not a matrix")}
  if (!is.matrix(Y)) { stop("Input X is not a matrix")}
  if (is_distance & !(is.matrix(X) & (nrow(X)==ncol(X)))) { stop("Input is not an accepted distance matrix") }
  if (!(is.logical(pca_center) && is.logical(pca_scale)) ) { stop("pca_center and pca_scale should be TRUE or FALSE")}
  if (!is.wholenumber(initial_dims) || initial_dims<=0) { stop("Incorrect initial dimensionality.")}

  tsne.args <- .check_tsne_params(nrow(X), dims=dims, perplexity=perplexity, theta=theta, max_iter=max_iter, verbose=verbose,
                                   Y_init=Y_init, stop_lying_iter=stop_lying_iter, mom_switch_iter=mom_switch_iter,
                                   momentum=momentum, final_momentum=final_momentum, eta=eta, exaggeration_factor=exaggeration_factor,
                                  binding_force=binding_force,coupled_period=coupled_period,uncoupled_period=uncoupled_period)



  # Check for missing values
  X <- na.fail(X)
  Y <- na.fail(Y)

  cat("Performing cosine norm.")
  cos.norm = FALSE
  if(cos.norm) {
    l2 <- cosineNorm(Xshared, mode="l2norm")
    Xshared <- .apply_cosine_norm(Xshared, l2)

    l2 <- cosineNorm(Yshared, mode="l2norm")
    Yshared <- .apply_cosine_norm(Yshared, l2)

  }


  print("About to perform PCA...")


  # Apply PCA
  if (!is_distance) {
    if (pca) {
      if(verbose) cat("Performing PCA\n")
      if(partial_pca){
        if (!requireNamespace("irlba", quietly = TRUE)) {stop("Package \"irlba\" is required for partial PCA. Please install it.", call. = FALSE)}
        X <- irlba::prcomp_irlba(X, n = initial_dims, center = pca_center, scale = pca_scale)$x
        Y <- irlba::prcomp_irlba(Y, n = initial_dims, center = pca_center, scale = pca_scale)$x

      }else{
        if(verbose & min(dim(X))>2500) cat("Consider setting partial_pca=TRUE for large matrices\n")
        X <- prcomp(X, retx=TRUE, center = pca_center, scale. = pca_scale, rank. = initial_dims)$x
        Y <- prcomp(Y, retx=TRUE, center = pca_center, scale. = pca_scale, rank. = initial_dims)$x
      }
    }

    print("Finished PCA...")

    if (check_duplicates) {
      if (any(duplicated(X))) { stop("Remove duplicates before running TSNE.") }
      if (any(duplicated(Y))) { stop("Remove duplicates before running TSNE.") }

    }



    if(is.null(mat12) && is.null(mat21)) {

      print("Finding MNNs...")
      # TODO: Replace with Aaron's C++ version
      # TODO: Apply PCA to shared dimensions

      mnns <- batchelor::findMutualNN(Yshared,Xshared,nk1,nk2)
      print("Found MNNs.")


      Z1 <- sapply(seq_len(nrow(Y)),function(i) mat[(mat[,1]==i),2])  #nns of data2 in data1
      Z2 <- sapply(seq_len(nrow(X)),function(i) mat[(mat[,2]==i),1])  #nns of data1 in data2

      mat21<-do.call(rbind,lapply(Z1,
                                  function(x)
                                    if(length(x)==0) rep(0,nk1)
                                  else if(length(x)<nk1) c(x,rep(0,nk1-length(x)))
                                  else x)) # nns of data 2 in data1

      mat12<-do.call(rbind,lapply(Z2,
                                  function(x)
                                    if(length(x)==0) rep(0,nk2)
                                  else if(length(x)<nk2) c(x,rep(0,nk2-length(x)))
                                  else x))
      mat12[mat12==0]<-NaN
      mat21[mat21==0]<-NaN

    }


    if (normalize) {
      print("Normalising...")
      X <- normalize_input(X)
      Y <- normalize_input(Y)

    }
    X <- t(X) # transposing for rapid column-major access.
    Y <- t(Y) # transposing for rapid column-major access.
  } else {
    # Compute Squared distance if we are using exact TSNE
    if (theta==0.0) {
      X <- X^2
      Y <- Y^2
    }
  }

  print("About to call Rsatsne_cpp code")
  out <- do.call(Rsatsne_cpp, c(list(X1=X, X2=Y, mat12, mat21, nk1, nk2, distance_precomputed=is_distance, num_threads=num_threads), tsne.args))

}



find.mutual.nn <- function(data1, data2, k1, k2)
  # Finds mutal neighbors between data1 and data2.
{

  n1 <- nrow(data1)
  n2 <- nrow(data2)
  n.total <- n1 + n2

  W21 <- FNN::get.knnx(data2, query=data1, k=k1)
  W12 <- FNN::get.knnx(data1, query=data2, k=k2)

  inds12<-cbind( as.vector(rep(seq_len(n1), each=k1)),  as.vector(t(W21$nn.index)) )
  inds21<-cbind(as.vector(t(W12$nn.index)) , as.vector(rep(seq_len(n2), each=k2)) )

  A<-rbind(inds12,inds21)
  keeps=duplicated(A)  ##duplicated rows of A are mutaul nns
  A<-A[keeps,]

  A1<-A[,1]
  A2<-A[,2]
  # Report cells that are MNNs.
  list(mnns=A)
}


