#include /include/common_functions.stan

data {
  // data
  // are individual observation or summary statistics used
  int is_ss;
  // sample sizes
  int<lower=0> N1;
  int<lower=0> N2;
  // individual observations
  vector[is_ss == 0 ? N1 : 0] x1;
  vector[is_ss == 0 ? N2 : 0] x2;
  // summary statistics
  vector[is_ss == 1 ? 2 : 0] mean_i;
  vector[is_ss == 1 ? 2 : 0] sd_i;

  // model type (1 if the parameter is included, 0 if not)
  int is_d;
  int is_r;
  int is_nu;

  // range of the parameters
  vector[is_d  == 1 ? 2 : 0] bounds_d;
  vector[is_r  == 1 ? 2 : 0] bounds_r;
  vector[is_nu == 1 ? 2 : 0] bounds_nu;
  array[is_d   == 1 ? 2 : 0] int bounds_type_d;
  array[is_r   == 1 ? 2 : 0] int bounds_type_r;
  array[is_nu == 1 ? 2 : 0] int bounds_type_nu;

  // prior distribution specification of the parameteres
  array[is_d   == 0 ? 1 : 0] real fixed_d;
  array[is_r   == 0 ? 1 : 0] real fixed_r;
  array[is_nu == 0 ? 1 : 0] real fixed_nu;
  vector[is_d  == 1 ? 3 : 0] prior_parameters_d;
  vector[is_r  == 1 ? 3 : 0] prior_parameters_r;
  vector[is_nu == 1 ? 3 : 0] prior_parameters_nu;
  int prior_type_d;
  int prior_type_r;
  int prior_type_nu;
}
parameters{
  real mu;
  real<lower = 0> sigma2;
  array[is_d] real<lower = coefs_lb(bounds_type_d,  bounds_d),  upper = coefs_ub(bounds_type_d,  bounds_d)>  delta;
  array[is_r] real<lower = coefs_lb(bounds_type_r,  bounds_r),  upper = coefs_ub(bounds_type_r,  bounds_r)>  rho;
  array[is_nu] real<lower = coefs_lb(bounds_type_nu, bounds_nu), upper = coefs_ub(bounds_type_nu, bounds_nu)> nu_p;
}
transformed parameters {
  real pooled_sigma;
  array[2] real sigma_i;
  array[2] real scale_i;
  array[2] real mu_i;
  real nu;


  // compute means and sigmas for each group
  if(is_nu == 1){
    nu = nu_p[1] + 2;
  }else{
    nu = fixed_nu[1];
  }
  if(is_r == 1){
    sigma_i[1]   = sqrt( 1 / (2 * 1/sigma2 * rho[1]       ) );
    sigma_i[2]   = sqrt( 1 / (2 * 1/sigma2 * (1 - rho[1]) ) );
    pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);
  }else{
    sigma_i[1]   = sqrt( 1 / (2 * 1/sigma2 * fixed_r[1]       ) );
    sigma_i[2]   = sqrt( 1 / (2 * 1/sigma2 * (1 - fixed_r[1]) ) );
    pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);
  }

  scale_i[1] = sigma_i[1] / sqrt(nu / (nu - 2.0));
  scale_i[2] = sigma_i[2] / sqrt(nu / (nu - 2.0));
  pooled_sigma = pool_sigma(sigma_i[1], sigma_i[2], N1, N2);

  if(is_d == 1){
    mu_i[1] = mu - 0.5 * delta[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * delta[1] * pooled_sigma;
  }else{
    mu_i[1] = mu - 0.5 * fixed_d[1] * pooled_sigma;
    mu_i[2] = mu + 0.5 * fixed_d[1] * pooled_sigma;
  }
}
model {
  // default Jeffrey's priors for mu and sigma
  target += Jeffreys_mu_lpdf(mu);
  target += Jeffreys_sigma_lpdf(sigma2);

  // priors on d and r
  if(is_d == 1){
    target += set_prior(delta[1], prior_type_d, prior_parameters_d, bounds_type_d, bounds_d);
  }
  if(is_r == 1){
    target += set_prior(rho[1], prior_type_r, prior_parameters_r, bounds_type_r, bounds_r);
  }
  if(is_nu == 1){
    target += set_prior(nu_p[1], prior_type_nu, prior_parameters_nu, bounds_type_nu, bounds_nu);
  }

  // likelihood of the data
  if(is_ss == 0){
    target += student_t_lpdf(x1 | nu, mu_i[1], scale_i[1]);
    target += student_t_lpdf(x2 | nu, mu_i[2], scale_i[2]);
  }else{
    reject("Fitting models with t likelihood and summary statistics is not possible :(.");
  }
}
