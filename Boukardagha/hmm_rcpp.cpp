// hmm_rcpp.cpp — C++ HMM engine for Boukardagha (2026) replication
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

static const double LOG2PI = std::log(2.0 * datum::pi);

// [[Rcpp::export]]
arma::vec dmvnorm_log_cpp(const arma::mat& X, const arma::vec& mu,
                           const arma::mat& Sigma) {
  int T = X.n_rows, d = X.n_cols;
  arma::mat L;
  if (!chol(L, Sigma, "lower"))
    chol(L, Sigma + 1e-6 * eye(d, d), "lower");
  double log_det = 2.0 * accu(log(L.diag()));
  double constant = -0.5 * (d * LOG2PI + log_det);
  arma::vec out(T);
  for (int t = 0; t < T; t++) {
    arma::vec z = X.row(t).t() - mu;
    arma::vec solved = solve(trimatl(L), z);
    out(t) = constant - 0.5 * dot(solved, solved);
  }
  return out;
}

static inline double log_sum_exp(const arma::vec& x) {
  double mx = x.max();
  if (!std::isfinite(mx)) return -datum::inf;
  return mx + std::log(accu(exp(x - mx)));
}

// [[Rcpp::export]]
List forward_backward_cpp(const arma::mat& B, const arma::mat& log_A,
                           const arma::vec& log_pi0) {
  int T = B.n_rows, K = B.n_cols;
  arma::mat log_alpha(T, K, fill::value(-datum::inf));
  arma::mat log_beta(T, K, fill::value(-datum::inf));
  arma::mat gamma(T, K);

  log_alpha.row(0) = (log_pi0 + B.row(0).t()).t();
  arma::vec work(K);
  for (int t = 1; t < T; t++)
    for (int k = 0; k < K; k++) {
      for (int j = 0; j < K; j++) work(j) = log_alpha(t-1, j) + log_A(j, k);
      log_alpha(t, k) = log_sum_exp(work) + B(t, k);
    }

  log_beta.row(T-1).fill(0.0);
  for (int t = T-2; t >= 0; t--)
    for (int j = 0; j < K; j++) {
      for (int k = 0; k < K; k++) work(k) = log_A(j, k) + B(t+1, k) + log_beta(t+1, k);
      log_beta(t, j) = log_sum_exp(work);
    }

  for (int t = 0; t < T; t++) {
    arma::vec lg = log_alpha.row(t).t() + log_beta.row(t).t();
    gamma.row(t) = exp(lg - log_sum_exp(lg)).t();
    for (int k = 0; k < K; k++)
      if (!std::isfinite(gamma(t, k)) || gamma(t, k) < 1e-300) gamma(t, k) = 1e-300;
    gamma.row(t) /= accu(gamma.row(t));
  }

  double ll = log_sum_exp(log_alpha.row(T-1).t());
  arma::mat A_num(K, K, fill::value(1e-10));
  for (int j = 0; j < K; j++)
    for (int k = 0; k < K; k++) {
      arma::vec xi(T-1);
      for (int t = 0; t < T-1; t++)
        xi(t) = log_alpha(t, j) + log_A(j, k) + B(t+1, k) + log_beta(t+1, k);
      double v = log_sum_exp(xi);
      if (std::isfinite(v)) A_num(j, k) = std::exp(v);
    }

  return List::create(Named("log_alpha")=log_alpha, Named("log_beta")=log_beta,
                      Named("gamma")=gamma, Named("ll")=ll, Named("A_num")=A_num);
}

// [[Rcpp::export]]
List fit_hmm_cpp(const arma::mat& X, List mu_init, List Sigma_init,
                  const arma::mat& A_init, const arma::vec& pi0_init,
                  int max_iter, double tol, double cov_reg) {
  int T = X.n_rows, d = X.n_cols, K = mu_init.size();
  std::vector<arma::vec> mu(K);
  std::vector<arma::mat> Sigma(K);
  for (int k = 0; k < K; k++) {
    mu[k] = as<arma::vec>(mu_init[k]);
    Sigma[k] = as<arma::mat>(Sigma_init[k]);
  }
  arma::mat A = A_init; arma::vec pi0 = pi0_init;
  arma::mat log_A = log(A + 1e-300);
  double prev_ll = -datum::inf, ll = -datum::inf;
  int iter = 0;
  arma::mat B(T, K), gamma(T, K);

  for (iter = 0; iter < max_iter; iter++) {
    for (int k = 0; k < K; k++) B.col(k) = dmvnorm_log_cpp(X, mu[k], Sigma[k]);
    B.elem(find_nonfinite(B)).fill(-1e10);
    List fb = forward_backward_cpp(B, log_A, log(pi0 + 1e-300));
    gamma = as<arma::mat>(fb["gamma"]);
    ll = as<double>(fb["ll"]);
    arma::mat A_num = as<arma::mat>(fb["A_num"]);
    if (!std::isfinite(ll)) break;
    if (std::abs(ll - prev_ll) < tol && iter > 1) { iter++; break; }
    prev_ll = ll;

    pi0 = clamp(gamma.row(0).t(), 1e-6, 1.0); pi0 /= accu(pi0);
    A = clamp(A_num, 1e-10, datum::inf);
    for (int j = 0; j < K; j++) A.row(j) /= accu(A.row(j));
    log_A = log(A + 1e-300);

    for (int k = 0; k < K; k++) {
      double ws = accu(gamma.col(k));
      if (!std::isfinite(ws) || ws < 1e-6) continue;
      arma::vec nm(d, fill::zeros);
      for (int t = 0; t < T; t++) nm += gamma(t, k) * X.row(t).t();
      nm /= ws;
      if (nm.has_inf() || nm.has_nan()) nm = mean(X, 0).t();
      mu[k] = nm;
      arma::mat ns(d, d, fill::zeros);
      for (int t = 0; t < T; t++) {
        arma::vec df = X.row(t).t() - mu[k];
        ns += gamma(t, k) * (df * df.t());
      }
      Sigma[k] = 0.5 * (ns / ws + (ns / ws).t()) + cov_reg * eye(d, d);
    }
  }

  // Filtered probabilities
  for (int k = 0; k < K; k++) B.col(k) = dmvnorm_log_cpp(X, mu[k], Sigma[k]);
  B.elem(find_nonfinite(B)).fill(-1e10);
  arma::mat la = as<arma::mat>(forward_backward_cpp(B, log_A, log(pi0+1e-300))["log_alpha"]);
  arma::mat filtered(T, K);
  for (int t = 0; t < T; t++) {
    arma::vec a = la.row(t).t();
    filtered.row(t) = exp(a - log_sum_exp(a)).t();
    for (int k = 0; k < K; k++)
      if (!std::isfinite(filtered(t,k))) filtered(t,k) = 1.0/K;
    filtered.row(t) /= accu(filtered.row(t));
  }

  List mo(K), so(K);
  for (int k = 0; k < K; k++) { mo[k]=mu[k]; so[k]=Sigma[k]; }
  return List::create(Named("pi0")=pi0, Named("A")=A, Named("mu")=mo,
                      Named("Sigma")=so, Named("K")=K, Named("d")=d,
                      Named("log_lik")=ll, Named("gamma")=gamma,
                      Named("filtered")=filtered, Named("n_iter")=iter);
}

// ─── Score a sequence under fitted HMM params (forward pass only) ───────────
// Used for Python-style model selection: score(full) - score(train)
// [[Rcpp::export]]
double score_sequence_cpp(const arma::mat& X, List mu_list, List Sigma_list,
                           const arma::mat& A, const arma::vec& pi0) {
  int T = X.n_rows, K = mu_list.size();
  arma::mat B(T, K);
  for (int k = 0; k < K; k++)
    B.col(k) = dmvnorm_log_cpp(X, as<arma::vec>(mu_list[k]), as<arma::mat>(Sigma_list[k]));
  B.elem(find_nonfinite(B)).fill(-1e10);

  arma::mat log_A = log(A + 1e-300);
  arma::vec log_pi0 = log(pi0 + 1e-300);
  arma::mat log_alpha(T, K, fill::value(-datum::inf));
  log_alpha.row(0) = (log_pi0 + B.row(0).t()).t();
  arma::vec work(K);
  for (int t = 1; t < T; t++)
    for (int k = 0; k < K; k++) {
      for (int j = 0; j < K; j++) work(j) = log_alpha(t-1, j) + log_A(j, k);
      log_alpha(t, k) = log_sum_exp(work) + B(t, k);
    }
  return log_sum_exp(log_alpha.row(T-1).t());
}

// [[Rcpp::export]]
arma::vec hmm_filter_step_cpp(const arma::vec& alpha_prev, const arma::vec& x,
                                List mu_list, List Sigma_list, const arma::mat& A) {
  int K = alpha_prev.n_elem;
  arma::vec pred = A.t() * alpha_prev;
  pred = clamp(pred, 1e-300, 1.0);
  arma::mat xm(1, x.n_elem); xm.row(0) = x.t();
  arma::vec le(K);
  for (int k = 0; k < K; k++)
    le(k) = dmvnorm_log_cpp(xm, as<arma::vec>(mu_list[k]), as<arma::mat>(Sigma_list[k]))(0);
  arma::vec la = log(pred) + le;
  double mx = la.max();
  if (!std::isfinite(mx)) return arma::vec(K, fill::value(1.0/K));
  arma::vec a = exp(la - mx); a /= accu(a);
  for (int k = 0; k < K; k++) if (!std::isfinite(a(k))) a(k) = 1.0/K;
  return a;
}


// ─── MVO with exact L1 via proximal gradient (ISTA) ─────────────────────────
// max mu'w - gamma*w'Sigma*w - tau*||w - w_prev||_1
// s.t. sum(w)=1, 0<=w<=w_max
// [[Rcpp::export]]
arma::vec solve_mvo_cpp(const arma::vec& mu, const arma::mat& Sigma,
                         const arma::vec& w_prev,
                         double gamma, double tau, double w_max,
                         int max_iter = 2000, double tol = 1e-10) {
  int N = mu.n_elem;

  // PSD projection
  arma::mat S = 0.5 * (Sigma + Sigma.t());
  arma::vec eigval; arma::mat eigvec;
  eig_sym(eigval, eigvec, S);
  eigval = clamp(eigval, 1e-10, datum::inf);
  S = eigvec * diagmat(eigval) * eigvec.t();

  // Step size from Lipschitz constant
  double L = 2.0 * gamma * eigval.max();
  if (L < 1e-12) L = 1.0;
  double step = 0.9 / L;

  arma::vec w = w_prev;

  for (int iter = 0; iter < max_iter; iter++) {
    // Gradient ascent on smooth part: mu - 2*gamma*Sigma*w
    arma::vec grad = mu - 2.0 * gamma * S * w;
    arma::vec z = w + step * grad;

    // Proximal for tau*||. - w_prev||_1: soft-thresholding around w_prev
    arma::vec diff = z - w_prev;
    arma::vec shrunk = sign(diff) % max(abs(diff) - step * tau, arma::zeros(N));
    arma::vec w_prox = w_prev + shrunk;

    // Project onto simplex with box constraints: sum=1, 0<=w<=w_max
    w_prox = clamp(w_prox, 0.0, w_max);
    // Duchi et al. simplex projection
    arma::vec u = sort(w_prox, "descend");
    arma::vec cssv = cumsum(u) - 1.0;
    int rho = 0;
    for (int j = N - 1; j >= 0; j--) {
      if (u(j) - cssv(j) / (j + 1.0) > 0) { rho = j; break; }
    }
    double theta = cssv(rho) / (rho + 1.0);
    arma::vec w_new = clamp(w_prox - theta, 0.0, w_max);
    // Re-normalize after box clamp (may have violated sum=1)
    double s = accu(w_new);
    if (s < 1e-10) {
      w_new.fill(1.0 / N);
    } else {
      w_new /= s;
    }

    if (arma::max(abs(w_new - w)) < tol) { w = w_new; break; }
    w = w_new;
  }

  return w;
}
