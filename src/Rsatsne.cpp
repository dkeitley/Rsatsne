#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

#include <Rcpp.h>
#include "tsne.h"
using namespace Rcpp;

Rcpp::List save_results(int N, int no_dims, const std::vector<double>& Y, const std::vector<double>& costs, const std::vector<double>& itercosts,
        double theta, double perplexity, int D, int stop_lying_iter, int mom_switch_iter,
        double momentum, double final_momentum, double eta, double exaggeration_factor);




arma::vec sample_mnns(NumericMatrix mnn_mat, unsigned int N) {
  Rcpp::NumericVector match (N);

  for(unsigned int i =0; i < N; i++) {
    arma::vec r = mnn_mat(i, Rcpp::_ );
    arma::uvec not_na = arma::find_finite(r);

    if(not_na.is_empty()) {
      match[i] = NA_REAL;

    } else {
      arma::uvec samp = Rcpp::RcppArmadillo::sample(not_na,1,false);
      match[i] = samp(0);
    }

  }

  arma::vec match_out = arma::vec(match);

  return(match_out);
}



arma::vec find_nearest_mnns(arma::mat ymat1, arma::mat ymat2, NumericMatrix mat21,
                       unsigned int N1, unsigned int N2,int no_dims, unsigned int nk1, unsigned int nk2) {


  arma::mat mnn_mat = Rcpp::as<arma::mat>(mat21);
  arma::mat mnn_emb_pos = arma::mat(N2*nk1,2).fill(NA_REAL);

  int start=0;
  for(int i=0; i<nk1; i++) {
    arma::uvec mnns_i = arma::conv_to<arma::uvec>::from(mnn_mat.col(i));

    arma::uvec not_na = find(mnns_i != 0);
    mnns_i = mnns_i.elem(not_na);

    arma::uvec inds = not_na + start;
    mnn_emb_pos.rows(inds) = ymat1.rows(mnns_i-1); //-1 to fix for C++ indexing

    start = start+N2;
  }


  arma::mat yd2 = arma::repmat(ymat2,nk1,1); //dataset 2 embedding positions


  // Calculate embedding distance between MNNs
  arma::mat ydiff = yd2-mnn_emb_pos;
  arma::mat ydiff2 = arma::square(ydiff);
  arma::vec lnorm = arma::sum(ydiff2,1);
  arma::mat tsdist12(lnorm);
  tsdist12.reshape(N2,nk1);

  // convert NAs to Inf
  arma::uvec na_vals = arma::find_nonfinite(tsdist12);
  tsdist12.elem(na_vals).fill(arma::datum::inf);
  arma::mat na_mat = arma::mat(N2,1).fill(arma::datum::inf);

  // Insert column of Infs to determine NA rows
  tsdist12.insert_cols(0,na_mat);

  arma::uvec closest_mnns = index_min(tsdist12,1);
  arma::uvec na_mnns = arma::find(closest_mnns==0);

  closest_mnns = closest_mnns - 1; // switch back to C++ indices

  arma::vec match = arma::conv_to<arma::vec>::from(closest_mnns);
  match.elem(na_mnns).fill(arma::datum::nan);


  return(match);

}



double * compute_mnn_gradient(double* Y1, double* Y2, NumericMatrix mat21, unsigned int nk1,
                                unsigned int N1, unsigned int N2, int no_dims, double mnn_weight) {

  // Load Y matrices by column
  arma::mat ymat1 = arma::mat(N1,no_dims);

  for(unsigned int i=0; i<(N1*no_dims);i+=no_dims) {
    for(unsigned int j=0; j<no_dims; j++) {
      ymat1(i/no_dims,j) = Y1[i+j];
    }
  }

  arma::mat ymat2 = arma::mat(N2,no_dims);

  for(unsigned int i=0; i<(N2*no_dims);i+=no_dims) {
    for(unsigned int j=0; j<no_dims; j++) {
      ymat2(i/no_dims,j) = Y2[i+j];
    }
  }

  arma::mat mnn_mat = Rcpp::as<arma::mat>(mat21);
  arma::mat mnn_emb_pos = arma::mat(N2*nk1,2).fill(NA_REAL);

  int start=0;
  for(int i=0; i<nk1; i++) {
    arma::uvec mnns_i = arma::conv_to<arma::uvec>::from(mnn_mat.col(i));

    arma::uvec not_na = find(mnns_i != 0);
    mnns_i = mnns_i.elem(not_na);

    arma::uvec inds = not_na + start;
    mnn_emb_pos.rows(inds) = ymat1.rows(mnns_i-1); //-1 to fix for C++ indexing

    start = start+N2;
  }


  arma::mat yd2 = arma::repmat(ymat2,nk1,1);
  arma::mat ydiff = yd2-mnn_emb_pos;
  arma::vec lnorm = arma::sum(ydiff,1);
  lnorm = mnn_weight*lnorm;

  return(lnorm.memptr());

}







double* calc_new_momentum(arma::mat ymat1,arma::mat ymat2, arma::vec match21,
                         NumericMatrix mat12, NumericMatrix mat21,
                         unsigned int N1, unsigned int N2, int no_dims,double binding_force) {


  // Get indices of chosen mnn in embedding matrix
  arma::uvec ydata2_inds = arma::uvec(N1);
  for(unsigned int i=0; i<N1; i++){
    ydata2_inds[i] = mat12(i,match21(i));
  }

  arma::uvec not_na2 = find(ydata2_inds != 0);

  arma::mat direct1 = arma::zeros(N1,no_dims);
  arma::mat direct2 = arma::zeros(N2,no_dims);

  direct1.rows(not_na2) = ymat2.rows(ydata2_inds.elem(not_na2)-1);
  arma::mat mom_direct_1 = (-ymat1 + direct1);
  mom_direct_1.rows(find(ydata2_inds == 0)).fill(0);


  // Return momentum as array in row order
  double* mom1 = (double*) malloc(N1 * no_dims * sizeof(double));

  for(unsigned int i=0; i<(N1*no_dims);i+=no_dims) {
    for(unsigned int j=0; j<no_dims; j++) {
      mom1[i+j] = binding_force*mom_direct_1(i/no_dims,j);
    }
  }

  return(mom1);

}




double* satsne_update(int iter, double* Y1,double*Y2, NumericMatrix mat12,NumericMatrix mat21,
                   unsigned int N1, unsigned int N2, int no_dims, int nk1, int nk2, bool sample, double binding_force) {


  // Load Y matrices by column
  arma::mat ymat1 = arma::mat(N1,no_dims);

  for(unsigned int i=0; i<(N1*no_dims);i+=no_dims) {
    for(unsigned int j=0; j<no_dims; j++) {
      ymat1(i/no_dims,j) = Y1[i+j];
    }
  }

  arma::mat ymat2 = arma::mat(N2,no_dims);

  for(unsigned int i=0; i<(N2*no_dims);i+=no_dims) {
    for(unsigned int j=0; j<no_dims; j++) {
      ymat2(i/no_dims,j) = Y2[i+j];
    }
  }

  arma::vec match12;

  // A) in early iterrations couple randomly with one of the sharedX MNNs
  if(sample) {
    match12 = sample_mnns(mat21,N2);

  // B) in late iterrations couple with the MNNs which is the closeset in embedding space
  } else {
     match12 = find_nearest_mnns(ymat1,ymat2,mat21,N1,N2,no_dims,nk1,nk2);
  }

  double* new_mom_2 = calc_new_momentum(ymat2,ymat1,match12,mat21,mat12,N2,N1,no_dims,binding_force);

  return(new_mom_2);

}



void rescale_embedding(double* Y, int N, int no_dims, int box_size) {

  //TODO: Generalise for no_dims
  double mins[2] = {Y[0],Y[1]};
  double maxs[2] = {Y[0],Y[1]};

  for(int j=0; j<no_dims; j++) {
    for(int i=0; i<(N*no_dims);i+=no_dims) {
      double cur_val = Y[i+j];
      if(cur_val < mins[j]) {mins[j] = cur_val;}
      if(cur_val > maxs[j]) {maxs[j] = cur_val;}
    }
  }


  for(int j=0; j<no_dims; j++) {
    for(int i=0; i<(N*no_dims);i+=no_dims) {
      Y[i+j] = box_size*((Y[i+j]-mins[j])/(maxs[j]-mins[j]));
    }
  }


}



void run_satsne(TSNE<2> tsne1, TSNE<2> tsne2, NumericMatrix mat12,NumericMatrix mat21,
                int nk1, int nk2, double* X1, int N1, int D1, double* Y1, double* X2, int N2,
                int D2, double * Y2, bool distance_precomputed, double* costs1,
                double* costs2, double* itercosts,int max_iter, double binding_force,
                int coupled_period, int uncoupled_period) {

  tsne1.initialise(X1, N1, D1, Y1, distance_precomputed, costs1, itercosts);
  tsne2.initialise(X2, N2, D2, Y2, distance_precomputed, costs2, itercosts);

  Rprintf("Initialised both tsnes.\n");

  int no_dims = 2;

  double* zero_mom1 = (double*) calloc(N1*no_dims, sizeof(double));

  //TODO: REMOVE THIS
  for(int i=0; i<(N1*no_dims); i++) {
    zero_mom1[i] = 0;
  }

  double* zero_mom2 = (double*) calloc(N2*no_dims, sizeof(double));
  for(int i=0; i<(N2*no_dims); i++) {
    zero_mom2[i] = 0;
  }


  int period = coupled_period + uncoupled_period;

  double* mom_update1;
  double* mom_update2;

  bool sample_mnns = false;

  max_iter= max_iter  - (max_iter % period ) + coupled_period + 10;

  Rprintf("Starting iterations...\n");

  for(int iter = 0; iter < max_iter; iter++) {

    if(iter<(max_iter/4)) {sample_mnns=true;}
    else {sample_mnns=false;}

    if(((iter+1) % period) < coupled_period) {
      mom_update2 = satsne_update(iter,Y1,Y2,mat12,mat21,N1,N2,2,nk1,nk2,sample_mnns,binding_force);
      mom_update1 = satsne_update(iter,Y2,Y1,mat21,mat12,N2,N1,2,nk1,nk2,sample_mnns,binding_force);

    } else {
      mom_update2 = zero_mom2;
      mom_update1 = zero_mom1;

     }

    double* mnn_grad2 = compute_mnn_gradient(Y1,Y2, mat21, nk1,N1,N2,no_dims, 1);
    double* mnn_grad1 = compute_mnn_gradient(Y2,Y1, mat12, nk2,N2,N1,no_dims, 1);



    tsne1.iterate(iter,N1,Y1,mom_update1,mnn_grad1,costs1,itercosts);
    tsne2.iterate(iter,N2,Y2,mom_update2,mnn_grad2,costs2,itercosts);


    int box_size = 10;
    //rescale_embedding(Y1,N1,no_dims,box_size);
    //rescale_embedding(Y2,N2,no_dims,box_size);

    //Rcpp::Rcout << std::endl << mom_update1[62] << std::endl;
    //Rcpp::Rcout << std::endl << mom_update2[132] << std::endl;


  }

}







// Function that runs the Barnes-Hut implementation of t-SNE
// [[Rcpp::export]]
Rcpp::List Rsatsne_cpp(NumericMatrix X1, NumericMatrix X2, NumericMatrix mat12,NumericMatrix mat21,
                       int nk1, int nk2, int no_dims, double perplexity,
                       double theta, bool verbose, int max_iter,
                     bool distance_precomputed, NumericMatrix Y_in, bool init,
                     int stop_lying_iter, int mom_switch_iter,
                     double momentum, double final_momentum, double binding_force,
                     int coupled_period, int uncoupled_period,
                     double eta, double exaggeration_factor, unsigned int num_threads) {

    size_t N1 = X1.ncol(), D1 = X1.nrow();
    size_t N2 = X2.ncol(), D2 = X2.nrow();

    double * data1=X1.begin();
    double * data2=X2.begin();

    if (verbose) {
      Rprintf("Read the %i x %i data matrix successfully!\n", N1, D1);
      Rprintf("Read the %i x %i data matrix successfully!\n", N2, D2);
    }
    std::vector<double> Y1(N1 * no_dims), costs1(N1), itercosts(static_cast<int>(std::ceil(max_iter/50.0)));
    std::vector<double> Y2(N2 * no_dims), costs2(N2);


    // Providing user-supplied solution.
    if (init) {
        for (size_t i = 0; i < Y1.size(); i++) Y1[i] = Y_in[i];
        for (size_t i = 0; i < Y2.size(); i++) Y2[i] = Y_in[i];
        if (verbose) Rprintf("Using user supplied starting positions\n");
    }


    // TODO: Add functionality for 1 or 3 dimensions output
    // Run tsne
    if (no_dims==1) {
      //SATSNE<1> satsne(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
        //      momentum, final_momentum, eta, exaggeration_factor, num_threads);
      //satsne.run(data, N1, D1, Y1.data(), distance_precomputed, costs.data(), itercosts.data());
    } else if (no_dims==2) {

      TSNE<2> tsne1(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
                    momentum, final_momentum, eta, exaggeration_factor, num_threads);

      TSNE<2> tsne2(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
                    momentum, final_momentum, eta, exaggeration_factor, num_threads);

      run_satsne(tsne1,tsne2,mat12,mat21,nk1,nk2, data1, N1, D1, Y1.data(), data2, N2,D2,Y2.data(),
                 distance_precomputed, costs1.data(),costs2.data(), itercosts.data(),max_iter,
                 binding_force,coupled_period, uncoupled_period);


    } else if (no_dims==3) {
      //TSNE<3> satsne(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
            //  momentum, final_momentum, eta, exaggeration_factor, num_threads);
      //satsne.run(data1, N1, D1, Y1.data(), distance_precomputed, costs1.data(), itercosts.data());
    } else {
      Rcpp::stop("Only 1, 2 or 3 dimensional output is suppported.\n");
    }

    return Rcpp::List::create(Rcpp::_["Y1"]=Rcpp::NumericMatrix(no_dims, N1, Y1.data()),
            Rcpp::_["Y2"]=Rcpp::NumericMatrix(no_dims, N2, Y2.data()),
            Rcpp::_["costs1"]=Rcpp::NumericVector(costs1.begin(), costs1.end()),
            Rcpp::_["costs2"]=Rcpp::NumericVector(costs2.begin(), costs2.end()),
            Rcpp::_["itercosts"]=Rcpp::NumericVector(itercosts.begin(), itercosts.end()));

}


/*
// Function that runs the Barnes-Hut implementation of t-SNE on nearest neighbor results.
// [[Rcpp::export]]
Rcpp::List Rtsne_nn_cpp(IntegerMatrix nn_dex, NumericMatrix nn_dist,
                     int no_dims, double perplexity,
                     double theta, bool verbose, int max_iter,
                     NumericMatrix Y_in, bool init,
                     int stop_lying_iter, int mom_switch_iter,
                     double momentum, double final_momentum,
                     double eta, double exaggeration_factor, unsigned int num_threads) {

    size_t N = nn_dex.ncol(), K=nn_dex.nrow(); // transposed - columns are points, rows are neighbors.
    if (verbose) Rprintf("Read the NN results for %i points successfully!\n", N);
    std::vector<double> Y(N * no_dims), costs(N), itercosts(static_cast<int>(std::ceil(max_iter/50.0)));

    // Providing user-supplied solution.
    if (init) {
        for (size_t i = 0; i < Y.size(); i++) Y[i] = Y_in[i];
        if (verbose) Rprintf("Using user supplied starting positions\n");
    }

    // Run tsne
    if (no_dims==1) {
      TSNE<1> tsne(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
              momentum, final_momentum, eta, exaggeration_factor, num_threads);
      tsne.run(nn_dex.begin(), nn_dist.begin(), N, K, Y.data(), costs.data(), itercosts.data());
    } else if (no_dims==2) {
      TSNE<2> tsne(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
              momentum, final_momentum, eta, exaggeration_factor, num_threads);
      tsne.run(nn_dex.begin(), nn_dist.begin(), N, K, Y.data(), costs.data(), itercosts.data());
    } else if (no_dims==3) {
      TSNE<3> tsne(perplexity, theta, verbose, max_iter, init, stop_lying_iter, mom_switch_iter,
              momentum, final_momentum, eta, exaggeration_factor, num_threads);
      tsne.run(nn_dex.begin(), nn_dist.begin(), N, K, Y.data(), costs.data(), itercosts.data());
    } else {
      Rcpp::stop("Only 1, 2 or 3 dimensional output is suppported.\n");
    }

    return Rcpp::List::create(Rcpp::_["Y"]=Rcpp::NumericMatrix(no_dims, N, Y.data()),
            Rcpp::_["costs"]=Rcpp::NumericVector(costs.begin(), costs.end()),
            Rcpp::_["itercosts"]=Rcpp::NumericVector(itercosts.begin(), itercosts.end()));
}
*/
