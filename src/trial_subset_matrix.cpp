#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// arma in; arma out
// [[Rcpp::export]]
mat submat_arma(arma::mat X, arma::colvec T) {
    mat y = X.rows(find(T == 1));
    return y;
}

// rcpp in; arma out
// [[Rcpp::export]]
mat submat_arma2(NumericMatrix X, NumericVector T) {
    mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
    colvec tIdx(T.begin(), T.size(), false); 
    mat y = Xmat.rows(find(tIdx == 1));
    return y;
}

// rcpp in; rcpp out
// [[Rcpp::export]]
NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) { 
    int n=X.nrow(), k=X.ncol();
    NumericMatrix out(sum(condition),k);
    for (int i = 0, j = 0; i < n; i++) {
        if(condition[i]) {
            out(j,_) = X(i,_);
            j = j+1;
        }
    }
    return(out);
}

/*** R
library("microbenchmark")

# simulate data
n=100
p=0.3
T=rbinom(n,1,p)
X=as.matrix(cbind(rnorm(n),rnorm(n)))

# compare output
identical(X[T==1,],submat_arma(X,T))
identical(X[T==1,],submat_arma2(X,T))
identical(X[T==1,],submat_rcpp(X,T))

# benchmark
microbenchmark(X[T==1,],submat_arma(X,T),submat_arma2(X,T),submat_rcpp(X,T),times=500)

# increase n
n=10000
p=0.3
T=rbinom(n,1,p)
X=as.matrix(cbind(rnorm(n),rnorm(n)))
# benchmark
microbenchmark(X[T==1,],submat_arma(X,T),submat_arma2(X,T),submat_rcpp(X,T),times=500)

*/