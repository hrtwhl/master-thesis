// =============================================================================
//  wasserstein_hmm.cpp
//
//  C++ kernels for the Wasserstein-HMM portfolio framework
//  (Boukardagha, "Explainable Regime-Aware Investing", 2026)
//
//  Exported to R:
//    - gaussian_hmm_em()           full Baum-Welch EM with full covariance
//    - gaussian_hmm_score()        log-likelihood of a sequence given a fitted model
//    - gaussian_hmm_predict_proba() smoothed posteriors gamma_t(k)
//    - wasserstein2_gaussian()     closed-form W2 between two Gaussians
//    - ledoit_wolf_cov()           Ledoit-Wolf shrinkage (scaled-identity target)
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// Helpers (file-local)
// -----------------------------------------------------------------------------

static inline mat sym_(const mat& A) { return 0.5 * (A + A.t()); }

// Symmetric PSD square root via eigendecomposition (robust to tiny negatives)
static mat sqrtm_psd_(const mat& A, double eps = 1e-12) {
  mat As = sym_(A);
  vec ev;
  mat V;
  if (!eig_sym(ev, V, As)) {
    return mat(A.n_rows, A.n_cols, fill::eye);
  }
  ev.transform([eps](double v) { return std::sqrt(std::max(v, eps)); });
  return V * diagmat(ev) * V.t();
}

template<typename T>
static inline double logsumexp_(const T& x) {
  double m = x.max();
  if (!std::isfinite(m)) return m;
  return m + std::log(arma::accu(exp(x - m)));
}

// Log multivariate Gaussian density at each row of X, for each component k.
// X: T x D, means: K x D, covars: D x D x K. Returns T x K.
static mat log_mvn_pdf_(const mat& X, const mat& means, const cube& covars,
                        double jitter = 1e-8) {
  uword T = X.n_rows, D = X.n_cols, K = means.n_rows;
  mat lp(T, K, fill::zeros);
  const double log2pi = std::log(2.0 * M_PI);
  
  for (uword k = 0; k < K; ++k) {
    mat Sk = sym_(covars.slice(k)) + jitter * eye(D, D);
    mat L;
    bool ok = chol(L, Sk, "lower");
    if (!ok) {
      Sk += 1e-4 * eye(D, D);
      ok = chol(L, Sk, "lower");
      if (!ok) { lp.col(k).fill(-datum::inf); continue; }
    }
    double log_det = 2.0 * arma::accu(log(L.diag()));
    rowvec mu = means.row(k);
    
    for (uword t = 0; t < T; ++t) {
      vec diff = (X.row(t) - mu).t();
      vec z = solve(trimatl(L), diff);          // L z = diff
      double q = arma::dot(z, z);               // (x-mu)' Sigma^{-1} (x-mu)
      lp(t, k) = -0.5 * (D * log2pi + log_det + q);
    }
  }
  return lp;
}

// K-means++ initialization
static mat kmeans_pp_(const mat& X, uword K, uint32_t seed) {
  uword T = X.n_rows, D = X.n_cols;
  mat C(K, D, fill::zeros);
  std::mt19937 gen(seed);
  std::uniform_int_distribution<uword> uni(0, T - 1);
  
  C.row(0) = X.row(uni(gen));
  vec d2(T, fill::value(datum::inf));
  
  for (uword k = 1; k < K; ++k) {
    for (uword t = 0; t < T; ++t) {
      double v = arma::accu(square(X.row(t) - C.row(k - 1)));
      if (v < d2(t)) d2(t) = v;
    }
    double total = arma::accu(d2);
    if (total <= 0) { C.row(k) = X.row(uni(gen)); continue; }
    
    std::uniform_real_distribution<double> ur(0.0, total);
    double r = ur(gen);
    double cs = 0.0;
    uword chosen = T - 1;
    for (uword t = 0; t < T; ++t) {
      cs += d2(t);
      if (cs >= r) { chosen = t; break; }
    }
    C.row(k) = X.row(chosen);
  }
  return C;
}

// Lloyd's k-means iterations starting from K-means++ centers
static void kmeans_lloyd_(const mat& X, uword K, mat& C, uvec& assign,
                          uword max_iter = 50) {
  uword T = X.n_rows, D = X.n_cols;
  assign.set_size(T);
  assign.zeros();
  
  for (uword it = 0; it < max_iter; ++it) {
    bool changed = false;
    for (uword t = 0; t < T; ++t) {
      double best = datum::inf;
      uword bi = 0;
      for (uword k = 0; k < K; ++k) {
        double d = arma::accu(square(X.row(t) - C.row(k)));
        if (d < best) { best = d; bi = k; }
      }
      if (assign(t) != bi) { changed = true; assign(t) = bi; }
    }
    mat Cnew(K, D, fill::zeros);
    vec cnt(K, fill::zeros);
    for (uword t = 0; t < T; ++t) {
      Cnew.row(assign(t)) += X.row(t);
      cnt(assign(t)) += 1.0;
    }
    for (uword k = 0; k < K; ++k) {
      if (cnt(k) > 0) Cnew.row(k) /= cnt(k);
      else Cnew.row(k) = C.row(k);     // keep old if empty
    }
    double shift = arma::norm(Cnew - C, "fro");
    C = Cnew;
    if (!changed || shift < 1e-7) break;
  }
}

// Stable forward pass in log domain.
// Returns log_alpha (TxK) and the joint log-likelihood log p(x_{1:T}).
static void forward_log_(const mat& log_b, const mat& log_A,
                         const vec& log_pi, mat& log_alpha, double& loglik) {
  uword T = log_b.n_rows, K = log_b.n_cols;
  log_alpha.set_size(T, K);
  log_alpha.row(0) = log_pi.t() + log_b.row(0);
  
  for (uword t = 1; t < T; ++t) {
    for (uword k = 0; k < K; ++k) {
      vec tmp = log_alpha.row(t - 1).t() + log_A.col(k);   // K-vec
      log_alpha(t, k) = logsumexp_(tmp) + log_b(t, k);
    }
  }
  loglik = logsumexp_(log_alpha.row(T - 1));
}

// Stable backward pass in log domain.
static void backward_log_(const mat& log_b, const mat& log_A, mat& log_beta) {
  uword T = log_b.n_rows, K = log_b.n_cols;
  log_beta.set_size(T, K);
  log_beta.row(T - 1).zeros();
  
  for (uword t = T - 1; t-- > 0; ) {
    for (uword k = 0; k < K; ++k) {
      vec tmp = log_A.row(k).t() + log_b.row(t + 1).t() + log_beta.row(t + 1).t();
      log_beta(t, k) = logsumexp_(tmp);
    }
  }
}

// -----------------------------------------------------------------------------
// 1) Gaussian HMM EM with full covariance.
//    Initialization: K-means++ -> Lloyd for emission means.
//    Transitions and start probabilities: uniform.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List gaussian_hmm_em(const arma::mat& X,
                     int K,
                     int max_iter = 100,
                     double tol = 1e-2,
                     int seed = 42,
                     double cov_reg = 1e-4) {
  
  uword T = X.n_rows, D = X.n_cols;
  if ((uword)K > T) stop("K must be <= number of observations");
  
  // --- Init via k-means
  mat C = kmeans_pp_(X, K, (uint32_t)seed);
  uvec assign;
  kmeans_lloyd_(X, K, C, assign, 50);
  
  mat means = C;                          // K x D
  cube covars(D, D, K, fill::zeros);
  vec pi_(K, fill::value(1.0 / K));
  mat A(K, K, fill::value(1.0 / K));      // K x K transitions
  
  // Initialize covariances from K-means clusters
  for (int k = 0; k < K; ++k) {
    uvec idx = find(assign == (uword)k);
    if (idx.n_elem > D + 1) {
      mat Xk = X.rows(idx);
      mat C0 = cov(Xk);                    // sample covariance
      covars.slice(k) = sym_(C0) + cov_reg * eye(D, D);
    } else {
      // Fallback: global covariance
      covars.slice(k) = sym_(cov(X)) + cov_reg * eye(D, D);
    }
  }
  
  mat log_A = log(A + 1e-300);
  vec log_pi = log(pi_ + 1e-300);
  
  double prev_ll = -datum::inf, ll = -datum::inf;
  int iter_done = 0;
  
  for (int it = 0; it < max_iter; ++it) {
    // --- E-step
    mat log_b = log_mvn_pdf_(X, means, covars);
    mat log_alpha, log_beta;
    forward_log_(log_b, log_A, log_pi, log_alpha, ll);
    backward_log_(log_b, log_A, log_beta);
    
    mat log_gamma = log_alpha + log_beta - ll;
    mat gamma_ = exp(log_gamma);
    
    // Vectorized xi_sum (sum_t over pairs (j,k))
    mat xi_sum(K, K, fill::zeros);
    for (uword t = 0; t < T - 1; ++t) {
      mat M = log_A;
      M.each_col() += log_alpha.row(t).t();
      M.each_row() += (log_b.row(t + 1) + log_beta.row(t + 1));
      xi_sum += exp(M - ll);
    }
    
    // --- M-step
    pi_ = gamma_.row(0).t() + 1e-300;
    pi_ /= arma::accu(pi_);
    log_pi = log(pi_);
    
    vec row_sums = sum(gamma_.rows(0, T - 2), 0).t();   // K-vec
    for (uword j = 0; j < (uword)K; ++j) {
      double rs = row_sums(j);
      if (rs < 1e-300) { A.row(j).fill(1.0 / K); continue; }
      A.row(j) = xi_sum.row(j) / rs;
      double s = arma::accu(A.row(j));
      if (s > 0) A.row(j) /= s;
    }
    log_A = log(A + 1e-300);
    
    // Vectorized emission means: mu_k = (gamma_k' X) / sum(gamma_k)
    vec gk_total = sum(gamma_, 0).t();
    for (uword k = 0; k < (uword)K; ++k) {
      double w = gk_total(k);
      if (w < 1e-300) continue;
      means.row(k) = (gamma_.col(k).t() * X) / w;
    }
    
    // Vectorized emission covariances:
    //   Sigma_k = ( (Xc .* gamma_k)' Xc ) / sum(gamma_k)
    for (uword k = 0; k < (uword)K; ++k) {
      double w = gk_total(k);
      if (w < 1e-300) {
        covars.slice(k) = sym_(cov(X)) + cov_reg * eye(D, D); continue;
      }
      mat Xc = X;
      Xc.each_row() -= means.row(k);                              // T x D centered
      mat Xc_w = Xc;
      Xc_w.each_col() %= gamma_.col(k);                           // weight each row
      mat S = (Xc_w.t() * Xc) / w;
      covars.slice(k) = sym_(S) + cov_reg * eye(D, D);
    }
    
    iter_done = it + 1;
    if (std::abs(ll - prev_ll) < tol) break;
    prev_ll = ll;
  }
  
  return List::create(
    _["means"]      = means,
    _["covars"]     = covars,
    _["transmat"]   = A,
    _["startprob"]  = pi_,
    _["log_lik"]    = ll,
    _["n_iter"]     = iter_done
  );
}

// -----------------------------------------------------------------------------
// 2) Log-likelihood of a NEW sequence given a fitted model.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double gaussian_hmm_score(const arma::mat& X,
                          const arma::mat& means,
                          const arma::cube& covars,
                          const arma::mat& transmat,
                          const arma::vec& startprob) {
  mat log_b = log_mvn_pdf_(X, means, covars);
  mat log_A = log(transmat + 1e-300);
  vec log_pi = log(startprob + 1e-300);
  mat log_alpha; double ll;
  forward_log_(log_b, log_A, log_pi, log_alpha, ll);
  return ll;
}

// -----------------------------------------------------------------------------
// 3) Smoothed posteriors gamma_t(k). Mirrors hmmlearn's predict_proba.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat gaussian_hmm_predict_proba(const arma::mat& X,
                                     const arma::mat& means,
                                     const arma::cube& covars,
                                     const arma::mat& transmat,
                                     const arma::vec& startprob) {
  mat log_b = log_mvn_pdf_(X, means, covars);
  mat log_A = log(transmat + 1e-300);
  vec log_pi = log(startprob + 1e-300);
  mat log_alpha, log_beta; double ll;
  forward_log_(log_b, log_A, log_pi, log_alpha, ll);
  backward_log_(log_b, log_A, log_beta);
  mat log_gamma = log_alpha + log_beta - ll;
  return exp(log_gamma);
}

// -----------------------------------------------------------------------------
// 4) Wasserstein-2 distance between two Gaussians (closed form).
//    W2^2 = ||mu1 - mu2||^2 + Tr( S1 + S2 - 2 ( S2^{1/2} S1 S2^{1/2} )^{1/2} )
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double wasserstein2_gaussian(const arma::vec& mu1, const arma::mat& S1,
                             const arma::vec& mu2, const arma::mat& S2) {
  mat S1s = sym_(S1), S2s = sym_(S2);
  double dm2 = arma::accu(square(mu1 - mu2));
  mat S2sqrt = sqrtm_psd_(S2s);
  mat M = S2sqrt * S1s * S2sqrt;
  mat Msqrt = sqrtm_psd_(M);
  double tr = arma::trace(S1s + S2s - 2.0 * Msqrt);
  if (tr < 0) tr = 0.0;
  return std::sqrt(dm2 + tr);
}

// -----------------------------------------------------------------------------
// 5) Ledoit-Wolf shrinkage with scaled-identity target (scikit-learn parity).
//    Returns the shrunk covariance matrix.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat ledoit_wolf_cov(const arma::mat& X) {
  uword n = X.n_rows, p = X.n_cols;
  if (n < 2) return eye(p, p);
  
  rowvec mu = mean(X, 0);
  mat Xc = X;
  Xc.each_row() -= mu;
  mat S = (Xc.t() * Xc) / (double)n;               // biased sample cov (sklearn)
  
  double mu_t = arma::trace(S) / (double)p;
  mat T_target = mu_t * eye(p, p);
  
  // delta^2 = ||S - mu*I||_F^2 / p
  double delta_ = arma::accu(square(S - T_target)) / (double)p;
  
  // beta^2 estimator (Eq. (10) in Chen et al. simplification):
  //   beta = ( sum_t ||x_t x_t' - S||_F^2 / n^2 ) / p
  //   shortcut: ( (1/n) * (X^2)' (X^2) - S^2 ).sum() / p
  mat Xc2 = square(Xc);
  mat E = (Xc2.t() * Xc2) / (double)n - square(S);
  double beta_ = arma::accu(E) / (double)p;
  
  if (beta_ > delta_) beta_ = delta_;
  if (beta_ < 0)      beta_ = 0;
  
  double shrinkage = (delta_ > 0) ? (beta_ / delta_) : 0.0;
  if (shrinkage > 1) shrinkage = 1;
  if (shrinkage < 0) shrinkage = 0;
  
  return (1.0 - shrinkage) * S + shrinkage * T_target;
}