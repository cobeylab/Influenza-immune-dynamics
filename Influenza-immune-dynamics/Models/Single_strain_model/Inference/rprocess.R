init <- Csnippet("
                 int N_STRAINS = n_strains;
                 int i,j;
                 double *I_vec = &I_s1; // initial infection status
                 double *h_vec = &h_s1; // initial titer 
                 double *h_last_infection_vec = &h_last_infection_s1; // initial titer at last infection 
                 double *t_last_infection_vec = &t_last_infection_s1; 
                 double *t_clear_vec = &t_clear_s1;
                 double *h_t0_vec = &h_t0_s1;
                 double *p_vec = &p_epi_1;
                 double *epi_vec = &epi_times_1;
                 double mean_gam = exp(log_mean_gam);
                 double var_gam = exp(log_var_gam);
                 int n_age_cat = n_age_categories;
                 double k;
                 double log_mean_boost;
                 double log_sd_boost;
                 double r = exp(log_r);
                 double T_peak = exp(log_T_peak);
                 double w = exp(log_w);
                 age = init_age;
                 double alpha;
                 double phi;
                 double infected = 0;
                 double d_general = exp(logit_d_general)/(1 + exp(logit_d_general));
                 double w_general; 
                
                if(h_baseline_individual_observed >= 10){
                  h_baseline_individual = runif(h_baseline_individual_observed, 2*h_baseline_individual_observed);
                }
                if(h_baseline_individual_observed < 10){
                  h_baseline_individual = runif(0,9);
                }

                d_acute = 0;
                if(age > age_thres){
                  k = exp(log_k_2);
                  alpha = alpha_2;
                  phi = phi_2;
                  w_general = exp(log_w_general_2);
                  log_mean_boost = log_mean_boost_2;
                  log_sd_boost = log_sd_boost_2;
                 }
                 if(age <= age_thres){
                  k = exp(log_k_1);
                  alpha = alpha_1;
                  phi = phi_1;
                  w_general = exp(log_w_general_1);
                  log_mean_boost = log_mean_boost_1;
                  log_sd_boost = log_sd_boost_1;
                 }
                 
            for(j = 0; j < N_STRAINS; j++){
                double size = 1;
                double p_temp;
                double p_sum = 0; 
                i = 0;
                while(size > 0){
                      if(i == 0){
                        p_temp = p_vec[i];
                      }
                      if(i > 0){
                        p_temp = p_vec[i]/(1 - p_sum);
                      }
                      double this_size = rbinom(size,p_temp);
                      size -= this_size;
                      p_sum += p_vec[i]; 
                      i++;
                }
                  
                t_last_infection_vec[j] = epi_vec[i];
                
                if(t_last_infection_vec[j] >= -5){
                   infected  = 1; 
                }
                if(infected == 1){
                    I_vec[j] = 1;
                    t_last_infection_vec[j] = t;
                    h_vec[j] = h_t0_vec[j];
                    h_last_infection_vec[j] = h_vec[j];
                    t_clear_vec[j] = t + rgamma((mean_gam * mean_gam)/var_gam, var_gam/mean_gam); 
                 }
                 if(infected == 0){
                    I_vec[j] = 0;
                    h_last_infection_vec[j] = h_baseline_individual;
                    if(variable_boosting == 1){
                        d_acute = exp(rnorm(log_mean_boost, exp(log_sd_boost)));
                    }
                    if(variable_boosting == 0){
                        d_acute = exp(log_mean_boost);
                    }
                    if(include_k == 1){
                      d_acute = d_acute*pow(h_last_infection_vec[j],-k);
                    }
                    double f_rise = h_baseline_individual*d_acute*(1-exp(-r*(T_peak)));
                    double h_peak = h_baseline_individual + f_rise;
                    double h_baseline = h_baseline_individual;
                    double f_wane =  (h_peak - h_baseline) * exp(-w*(t-(T_peak + t_last_infection_vec[j])));
                    h_vec[j]= h_baseline + f_wane;
                    t_clear_vec[j] = t_last_infection_vec[j] + rgamma((mean_gam * mean_gam)/var_gam, var_gam/mean_gam); 
                 }
                    
                    q1 =  1/(1 + exp(phi * (log(h_vec[j]) - alpha)));
                 
                    if(include_general_immunity == 1){
                        q2 = (1-(1-d_general)*exp(-w_general*((t - t_last_infection_vec[j]))));
                    }
                 
                    if(include_general_immunity == 0){
                      q2 = 1; 
                    }
                 
                    if(q1 > q2){
                      q = q2;
                    }
                    if(q1 <= q2){
                      q = q1;
                    }
                 }
                 ")

rprocess <- Csnippet("// set pointers for arrays
                     double *I_vec = &I_s1;
                     double *h_last_infection_vec = &h_last_infection_s1;
                     double *h_vec = &h_s1;
                     double *t_clear_vec = &t_clear_s1;
                     double *t_last_infection_vec = &t_last_infection_s1;
                     double *contact_group_age_vec = &age_contact_group_1;
                     double *beta_community_vec = &beta_community_1;

                     int N_STRAINS = n_strains;
                     int i,j;
                     
                     // transform params
                     double alpha; 
                     double phi; 
                     double infected;
                     double mean_gam = exp(log_mean_gam);
                     double var_gam = exp(log_var_gam);
                     double T_peak = exp(log_T_peak);
                     double r = exp(log_r);
                     double w = exp(log_w);
                     double beta_community;
                     double beta_scaled = exp(log_beta_scaled);
                     double k;
                     double log_mean_boost;
                     double log_sd_boost;
                     double s_group_1 = exp(log_imprinting_effect_group_1);
                     double s_group_2 = exp(log_imprinting_effect_group_2);
                     int n_age_cat = n_age_categories;
                     double d2;
                     double d_general = exp(logit_d_general)/(1 + exp(logit_d_general));
                     double w_general;

                    for(i = 0; i < (n_age_cat - 1); i++){
                      if(age >= contact_group_age_vec[i] & age < contact_group_age_vec[i + 1]){
                          beta_community = beta_community_vec[i];
                        }
                    }

                     if(age >= contact_group_age_vec[n_age_cat - 1]){
                        beta_community = beta_community_vec[n_age_cat - 1];
                     }
                      
                     for(j = 0; j < N_STRAINS; j++) { 
                        double I_temp = I_vec[j]; //track this strain
                        double h_temp;
                        double h_last_infection_temp = h_last_infection_vec[j];
                        double t_last_infection_temp = t_last_infection_vec[j];
                        double t_clear_temp = t_clear_vec[j];
                        double d_temp;

                       //choose titer params based on age and titer at last infection:
                       if(age > age_thres){
                        alpha = alpha_2;
                        phi = phi_2;
                        k = exp(log_k_2);
                        log_mean_boost = log_mean_boost_2;
                        log_sd_boost = log_sd_boost_2;
                        w_general = exp(log_w_general_2);
                        d2 = exp(logit_d2_2)/(1 + exp(logit_d2_2));
                       }
                       if(age <= age_thres){
                        alpha = alpha_1;
                        phi = phi_1;
                        k = exp(log_k_1);
                        log_mean_boost = log_mean_boost_1;
                        log_sd_boost = log_sd_boost_1;
                        w_general = exp(log_w_general_1);
                        d2 = exp(logit_d2_1)/(1 + exp(logit_d2_1));
                       }

                      //Update titer:
                       // If titer is currently rising
                       if(t <= (t_last_infection_temp + T_peak)){
                         double f_rise = h_last_infection_temp*d_acute*(1-exp(-r*(t-t_last_infection_temp)));
                         h_temp = h_last_infection_temp + f_rise;
                       } 
                       // If titer is currently waning 
                       if(t > (t_last_infection_temp + T_peak)){
                          double f_rise = h_last_infection_temp*d_acute*(1-exp(-r*(T_peak)));
                          double h_peak = h_last_infection_temp + f_rise;
                          double h_baseline = h_baseline_individual;
                          double f_wane =  (h_peak - h_baseline) * exp(-w*(t-(T_peak + t_last_infection_temp)));
                          h_temp = h_baseline + f_wane;
                        }
                                     
                     // determine FOI based on ab titer 
                    q1  =  1/(1 + exp(phi * (log(h_temp) - alpha)));
                    
                    if(include_general_immunity == 1){
                        q2 = (1-(1-d_general)*exp(-w_general*((t - t_last_infection_temp))));
                   }

                    if(include_general_immunity == 0){
                        q2 = 1; 
                    }

                    if(q1 > q2){
                      q = q2;
                    }
                    if(q1 <= q2){
                      q = q1;
                    }

                    // determine risk based on community FOI 

                    double lambda = q*(beta_community*beta_scaled*L);

                     // imprinting
                     if(include_imprinting == 1){
                       if(imprinting_group == 1){
                          lambda = lambda*(p_imprinted_group_1*s_group_1) + (1 - p_imprinted_group_1)*lambda;
                       }
                       if(imprinting_group == 2){
                          lambda = lambda*(p_imprinted_h3*s_group_2) + (1-p_imprinted_h3)*lambda;
                       }
                    }
                    
                    if(I_temp == 0){
                      double p_infection = (1-exp(-lambda*dt));
                      infected = (double)rbinom(1, p_infection);
                      if(infected == 1){
                        t_last_infection_temp = t;
                        h_last_infection_temp = h_temp;

                        //draw infection duration 
                        double duration = rgamma((mean_gam * mean_gam)/var_gam, var_gam/mean_gam); 
                        t_clear_temp = t + duration;
                        double h_baseline = h_baseline_individual;

                        if(variable_boosting == 1){
                          d_acute = exp(rnorm(log_mean_boost, exp(log_sd_boost)));
                        }
                        if(variable_boosting == 0){
                          d_acute = exp(log_mean_boost);
                        }
                        if(include_k == 1){
                          d_acute = d_acute*pow(h_last_infection_temp,-k);
                        }
                        h_baseline = h_baseline + h_baseline*d2*d_acute;
                        h_baseline_individual = h_baseline;
                      }
                     }
                     
                     if(I_temp == 1){
                       if(t >= t_clear_temp){
                         infected = 0;
                       }
                       if(t < t_clear_temp){
                         infected = 1;
                       }
                     }
                     
                     I_vec[j] = infected;
                     t_clear_vec[j] = t_clear_temp;
                     h_vec[j] = h_temp;
                     h_last_infection_vec[j] = h_last_infection_temp;
                     t_last_infection_vec[j] = t_last_infection_temp;
          
                   } // end strain loop
                     
                   // age this individual 
                   age = age + dt/365;
                  
                ")

## Measurement model ## --------------------------------------------------------------------------------
dmeasure <- Csnippet("
                     double sig = exp(log_sig);
                     double sig_2 = exp(log_sig_2);
                     double h_observed;
                     double h_latent;
                     double this_lik;
                     
                     // set pointers for arrays
                     double *I_vec = &I_s1;
                     double *h_vec = &h_s1;
                     double *t_clear_vec = &t_clear_s1;
                     double *h_obs_vec = &h_obs_s1;
            
                     int N_STRAINS = n_strains;
                     int j;
                     lik = 0;
                     for( j = 0; j < N_STRAINS; j ++){

                       if (ISNA(h_obs_vec[j])) {
                          lik += (give_log) ? 0 : 1;
                       } else {
                          if(log_transform_titers == 1){
                            h_latent = log2(h_vec[j]/10) + 2;
                          }
                          if(log_transform_titers == 0){
                            h_latent = h_vec[j];
                          }
                          if(log_transform_obs == 1){
                            h_observed = log2(h_obs_vec[j]/10) + 2;
                          }
                          if(log_transform_obs == 0){
                             h_observed = h_obs_vec[j];
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
                       }
                     
                     }
                     
                     lik = give_log ? lik : exp(lik);
                     ")

rmeasure <- Csnippet("
                     double sig = exp(log_sig);
                     double sig_2 = exp(log_sig_2);
                     // set pointers for arrays
                     double *I_vec = &I_s1;
                     double *h_vec = &h_s1;
                     double *h_obs_vec = &h_obs_s1;
                     int N_STRAINS = n_strains;
                     int j;
                     double h_latent;
                     
                     for( j = 0; j < N_STRAINS; j ++){
                        if(log_transform_titers == 1){
                          h_latent = log2(h_vec[j]/10) + 2;
                        }
                        if(log_transform_titers == 0){
                          h_latent = h_vec[j];
                        }
         
                          double rand;
                          if(h_latent < 2){
                            rand = rnorm(h_latent, sig_2);
                          }
                          if(h_latent >= 2){
                            rand = rnorm(h_latent,sig);
                          }
                          if(rand < 2){
                            rand = 1;
                          }
                          for(int z = 3; z < 11; z++){
                            double bound = (double)z;
                            if(rand > (bound-1)){
                              if(rand <= bound){
                                rand = bound - 1;
                              }
                            }
                          }
                          if(rand >= 11){
                            rand = 11;
                          }
                          h_obs_vec[j] = rand;
                     }
                     
                     ")
