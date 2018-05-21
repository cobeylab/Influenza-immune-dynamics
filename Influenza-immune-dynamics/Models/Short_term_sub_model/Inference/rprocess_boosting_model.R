init <- Csnippet("
                  h_init = runif(h_t0, 2*h_t0);
                  h = h_init;
                  age = init_age;
                 ")

rprocess <- Csnippet("
            
                     // transform params
                     double mean_gam = exp(log_mean_gam);
                     double var_gam = exp(log_var_gam);
                     double T_peak = exp(log_T_peak);
                     double r = exp(log_r);
                     double w = exp(log_w);
                     double k = exp(log_k);
                     double d_acute;
                     double t_inf = t_infection;

                     if(variable_boosting == 1){
                        d_acute = exp(rnorm(log_mean_boost, exp(log_sd_boost)));
                     }
                     if(variable_boosting == 0){
                        d_acute = exp(log_mean_boost);
                     }

                     if(t < t_inf){
                       h = h_init;
                     }
                     if(include_k == 1){
                      d_acute = d_acute*pow(h_init,-k);
                    }

                     if(t >=t_inf){
                       if(t <= (t_inf + T_peak)){
                         double f_rise = h_init*d_acute*(1-exp(-r*(t-t_inf)));
                         h = h_init + f_rise;
                        } 
                       if(t > (t_infection + T_peak)){
                          double f_rise = h_init*d_acute*(1-exp(-r*(T_peak)));
                          double h_peak = h_init + f_rise;
                          double f_wane =  (h_peak - h_init) * exp(-w*(t-(T_peak + t_inf)));
                          h = h_init + f_wane;
                        }
                     }
                   age = age + dt/365;
                  
                ")

## Measurement model ## --------------------------------------------------------------------------------
dmeasure <- Csnippet("
                     double sig = exp(log_sig);
                     double sig_2 = exp(log_sig_2);
                     double h_observed;
                     double h_latent;
                     double this_lik;
                     lik = 0;
                    
                     if (ISNA(h)) {
                          lik += (give_log) ? 0 : 1;
                       } else {
                          if(log_transform_titers == 1){
                            h_latent = log2(h/10) + 2;
                          }
                          if(log_transform_titers == 0){
                            h_latent = h;
                          }
                          if(log_transform_obs == 1){
                            h_observed = log2(h_obs/10) + 2;
                          }
                          if(log_transform_obs == 0){
                             h_observed = h_obs;
                          }
                 
                        if(h_observed < 2){
                            this_lik = log(pnorm(2,h_latent,sig_2,1,0));
                        }
                        if(h_observed >= 2){
                          if(h_observed < 11){
                            this_lik = log(pnorm(h_observed + 1,h_latent,sig,1,0) - pnorm(h_observed, h_latent,sig,1,0));
                          }
                          if(h_observed >= 11){
                            this_lik = log(1 - pnorm(11,h_latent,sig,1,0));
                          }
                        }
                     

                     lik += this_lik;
                     lik = give_log ? lik : exp(lik);
}
                     ")

rmeasure <- Csnippet("

                     double sig = exp(log_sig);
                     double sig_2 = exp(log_sig_2);
                     double h_latent;
                    
                        if(log_transform_titers == 1){
                          h_latent = log2(h/10) + 2;
                        }
                        if(log_transform_titers == 0){
                          h_latent = h;
                        }

                          double rand;
                          
                          if(h_latent < 2){
                            rand = rnorm(h_latent,sig_2);
                          }
                          if(h_latent >= 2){
                            rand = rnorm(h_latent,sig);
                          }
                      
                          if(rand <= 1){
                              rand = 1;
                            }
  
                            for(int k = 2; k < 11; k++){
                              double bound = (double)k;
                              if(rand > (bound-1)){
                                if(rand < bound){
                                  rand = bound - 1;
                                }
                              }
                            }
                            if(rand >= 11){
                              rand = 11;
                            }
                            h_obs = rand;

                     ")
