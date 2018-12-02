init <- Csnippet("
                 int N_MEMBERS = n_members;
                 int i,j;
                 int z;
                 double *I_vec = &I_m1; // initial infection status
                 double *h_vec = &h_m1; // initial titer 
                 double *h_last_infection_vec = &h_last_infection_m1; // initial titer at last infection 
                 double *t_last_infection_vec = &t_last_infection_m1; 
                 double *t_clear_vec = &t_clear_m1;
                 double *h_baseline_individual_vec = &h_baseline_individual_m1;
                 double *h_t0_vec = &h_t0_m1;
                 double *d_acute_vec = &d_acute_m1;
                
                 double *h_baseline_observed_vec = &h_baseline_observed_m1;
                 double *q1_vec = &q1_m1;
                 double *q2_vec = &q2_m1;
                 double *q_vec = &q_m1;
                 double *p_vec_1 = &p_epi_1_m1;
                 double *epi_vec_1 =  &epi_times_1_m1;
                 double *p_vec_2 = &p_epi_1_m2;
                 double *epi_vec_2 =  &epi_times_1_m2;
                 double *p_vec_3 = &p_epi_1_m3;
                 double *epi_vec_3 =  &epi_times_1_m3;
                 double *p_vec_4 = &p_epi_1_m4;
                 double *epi_vec_4 =  &epi_times_1_m4;
                 double *age_vec = &age_m1;
                 double *init_age_vec = &init_age_m1;

                 double mean_gam = exp(log_mean_gam);
                 double var_gam = exp(log_var_gam);
                 double r = exp(log_r);
                 double T_peak = exp(log_T_peak);
                 double w = exp(log_w);
                 double alpha;
                 double phi;
                 double infected = 0;
                 double d_general = exp(logit_d_general)/(1 + exp(logit_d_general));
                 double w_general; 
                 double k;
                 double log_mean_boost;
                 double log_sd_boost;
                 double weight;
          
                 for(j = 0; j < N_MEMBERS; j++){
                   age_vec[j] = init_age_vec[j];
                   if(age_vec[j] > age_thres){
                      alpha = alpha_adults;
                      k = exp(log_k_adults);
                      phi = phi_adults;
                      w_general = exp(log_w_general_adults);
                      log_mean_boost = log_mean_boost_adults;
                      log_sd_boost = log_sd_boost_adults;
                      weight = exp(logit_weight_adults)/(1 + exp(logit_weight_adults));
                   }
                   if(age_vec[j] <= age_thres){
                      alpha = alpha_kids;
                      k = exp(log_k_kids);
                      phi = phi_kids;
                      w_general = exp(log_w_general_kids);
                      log_mean_boost = log_mean_boost_kids;
                      log_sd_boost = log_sd_boost_kids;
                      weight = exp(logit_weight_kids)/(1 + exp(logit_weight_kids));
                   }

                   double q1 = q1_vec[j];
                   double q2 = q2_vec[j];
                   d_acute_vec[j] = 0;
                 
                  if(h_baseline_observed_vec[j] >= 10){
                    h_baseline_individual_vec[j] = runif(h_baseline_observed_vec[j], 2*h_baseline_observed_vec[j]);
                  }
                  if(h_baseline_observed_vec[j] < 10){
                   h_baseline_individual_vec[j] = runif(0,10);
                  }
                
                  double size = 1;
                  double p_temp;
                  double p_sum = 0; 
                  i = 0;

                  while(size > 0){
                      if(i == 0){
                        if(j == 0){
                          p_temp = p_vec_1[i];
                        }
                        if(j == 1){
                          p_temp = p_vec_2[i];
                          }
                        if(j == 2){
                          p_temp = p_vec_3[i];
                        }
                        if(j == 3){
                          p_temp = p_vec_4[i];
                        }
                      }
                      if(i > 0){
                        if(j == 0){
                          p_temp = p_vec_1[i]/(1 - p_sum);
                        }
                        if(j == 1){
                          p_temp = p_vec_2[i]/(1 - p_sum);
                        }
                        if(j == 2){
                          p_temp = p_vec_3[i]/(1 - p_sum);
                        }
                        if(j == 3){
                          p_temp = p_vec_4[i]/(1 - p_sum);
                        }
                      }

                      double this_size = rbinom(size,p_temp);
                      size -= this_size;
                      
                      if(j == 0){
                        p_sum += p_vec_1[i];
                      }
                      if(j == 1){
                        p_sum += p_vec_2[i]; 
                      }
                      if(j == 2){
                        p_sum += p_vec_3[i]; 
                      }
                      if(j == 3){
                        p_sum += p_vec_4[i]; 
                      }
                      i++;
                  }

                if( j == 0){
                  t_last_infection_vec[j] = epi_vec_1[i];
                 }
                 if( j == 1){
                 t_last_infection_vec[j] = epi_vec_2[i];
                 }
                if( j == 2){
                  t_last_infection_vec[j] = epi_vec_3[i];
                 }
                 if( j == 3){
                  t_last_infection_vec[j] = epi_vec_4[i];
                 }
                 
               
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
                    h_last_infection_vec[j] = h_baseline_individual_vec[j];
                    if(variable_boosting == 1){
                        d_acute_vec[j] = exp(rnorm(log_mean_boost, exp(log_sd_boost)));
                    }
                    if(variable_boosting == 0){
                        d_acute_vec[j] = exp(log_mean_boost);
                    }
                    if(include_k == 1){
                      d_acute_vec[j] = d_acute_vec[j]*pow(h_last_infection_vec[j],-k);
                    }
                    double f_rise = h_baseline_individual_vec[j]*d_acute_vec[j]*(1-exp(-r*(T_peak)));
                    double h_peak = h_baseline_individual_vec[j] + f_rise;
                    double h_baseline = h_baseline_individual_vec[j];
                    double f_wane =  (h_peak - h_baseline) * exp(-w*(t-(T_peak + t_last_infection_vec[j])));
                    h_vec[j]= h_baseline + f_wane;
                    t_clear_vec[j] = t_last_infection_vec[j] + rgamma((mean_gam * mean_gam)/var_gam, var_gam/mean_gam); 
                 }
                    
                    q1_vec[j] =  1/(1 + exp(phi * (log(h_vec[j]) - alpha)));
                 
                    if(include_general_immunity == 1){
                      q2_vec[j] = (1-(1-d_general)*exp(-w_general*((t - t_last_infection_vec[j]))));
                    }
                 
                    if(include_general_immunity == 0){
                      q2_vec[j] = 0; 
                    }
                   /* if(q1_vec[j] > q2_vec[j]){
                      q_vec[j] = q2_vec[j];
                    }
                    if(q1_vec[j] <= q2_vec[j]){
                      q_vec[j] = q1_vec[j];
                    }*/
                    q_vec[j] = q1_vec[j]*weight + q2_vec[j]*(1-weight);
          }
   
                 ")

rprocess <- Csnippet("// set pointers for arrays
                     double *I_vec = &I_m1;
                     double *h_last_infection_vec = &h_last_infection_m1;
                     double *h_vec = &h_m1;
                     double *t_clear_vec = &t_clear_m1;
                     double *t_last_infection_vec = &t_last_infection_m1;
                     double *h_baseline_individual_vec = &h_baseline_individual_m1;
                     double *d_acute_vec = &d_acute_m1;
                     double *contact_group_age_vec = &age_contact_group_1;
                     double *beta_community_vec = &beta_community_1;
                     double *L_vec = &L_m1;
                     double *q1_vec = &q1_m1;
                     double *q2_vec = &q2_m1;
                     double *q_vec = &q_m1;
                     double *age_vec = &age_m1;
                     double *p_imprinted_group_1_vec = &p_imprinted_group_1_m1;
                     double *p_imprinted_group_2_vec = &p_imprinted_h3_m1;

                     int N_MEMBERS = n_members;
                     int n_age_cat = n_age_categories;
                     int i,j;
                     int z;
                     double k;
                     double log_mean_boost;
                     double log_sd_boost;
                    
                     // transform params
                     double alpha; 
                     double phi; 
                     double infected;
                     double mean_gam = exp(log_mean_gam);
                     double var_gam = exp(log_var_gam);
                     double T_peak = exp(log_T_peak);
                     double omega_hh = exp(log_omega_hh);
                     double r = exp(log_r);
                     double w = exp(log_w);
                     double beta_community;
                     double s_group_1 = exp(log_imprinting_effect_group_1);
                     double s_group_2 = exp(log_imprinting_effect_group_2);
                     double d2;
                     double d_general = exp(logit_d_general)/(1 + exp(logit_d_general));
                     double beta_scaled = exp(log_beta_scaled);
                     double w_general;
                     double weight;

                
                     for(j = 0; j < N_MEMBERS; j++) { 
                        double lambda_community;
                        double lambda_household = 0;
                        double age = age_vec[j];
                        
                        for(i = 0; i < (n_age_cat - 1); i++){
                          if(age >= contact_group_age_vec[i] & age < contact_group_age_vec[i + 1]){
                            beta_community = beta_community_vec[i];
                          }
                        }
  
                       if(age >= contact_group_age_vec[n_age_cat - 1]){
                          beta_community = beta_community_vec[n_age_cat - 1];
                       }

                        double I_temp = I_vec[j];
                        double h_temp;
                        double h_last_infection_temp = h_last_infection_vec[j];
                        double t_last_infection_temp = t_last_infection_vec[j];
                        double h_baseline_individual = h_baseline_individual_vec[j];
                        double d_acute = d_acute_vec[j];
                        double t_clear_temp = t_clear_vec[j];
                        
                        double q1 = q1_vec[j];
                        double q2 = q2_vec[j];
                        double q = q_vec[j];
                        double L = L_vec[j];
                        
                       //choose titer params based on age and titer at last infection:
                       if(age > age_thres){
                          alpha = alpha_adults;
                          phi = phi_adults;
                          k = exp(log_k_adults);
                          w_general = exp(log_w_general_adults);
                          d2 = exp(logit_d2_adults)/(1 + exp(logit_d2_adults));
                          log_mean_boost = log_mean_boost_adults;
                          log_sd_boost = log_sd_boost_adults;
                          weight = exp(logit_weight_adults)/(1 + exp(logit_weight_adults));
                       }
                       if(age <= age_thres){
                          alpha = alpha_kids;
                          phi = phi_kids;
                          k = exp(log_k_kids);
                          w_general = exp(log_w_general_kids);
                          d2 = exp(logit_d2_kids)/(1 + exp(logit_d2_kids));
                          log_mean_boost = log_mean_boost_kids;
                          log_sd_boost = log_sd_boost_kids;
                          weight = exp(logit_weight_kids)/(1 + exp(logit_weight_kids));
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
                        q2 = 0; 
                    }

                    /*if(q1 > q2){
                      q = q2;
                    }
                    if(q1 <= q2){
                      q = q1;
                    }*/
                    q = q1*weight + q2*(1-weight);

                    // determine risk based on community FOI 
                    lambda_community = q*(beta_community*beta_scaled*L);
                    
                    for(z = 0; z < N_MEMBERS; z ++){
                      if(z != j){
                        if(I_vec[z] > 0){
                          lambda_household = q*omega_hh;
                        }
                      }
                    }
                    //Rprintf(\"t is %f\\n\",t);
                    //Rprintf(\"L is %f\\n\",L);
                    //Rprintf(\"member is %d\\n\",j);
                    //Rprintf(\"lambda_comm is %f\\n\",lambda_community);
                    //Rprintf(\"lambda_hh is %f\\n\",lambda_household);
                    
                    
                    double lambda = lambda_community + lambda_household;
                  
                    //Rprintf(\"lambda is %f\\n\",lambda);

                     // imprinting
                     if(include_imprinting == 1){
                       if(imprinting_group == 1){
                          lambda = lambda*(p_imprinted_group_1_vec[j]*s_group_1) + (1 - p_imprinted_group_1_vec[j])*lambda;
                       }
                       if(imprinting_group == 2){
                          lambda = lambda*(p_imprinted_group_2_vec[j]*s_group_2) + (1-p_imprinted_group_2_vec[j])*lambda;
                       }
                    }
                    
                    if(I_temp == 0){
                      double p_infection = (1-exp(-lambda*dt));
                      
                      //Rprintf(\"p_infection is %f\\n\",p_infection);
                      infected = (double)rbinom(1, p_infection);
                      if(infected == 1){
                        I_vec[j] = 1;
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
                          d_acute = d_acute*pow(h_last_infection_vec[j],-k);
                        }

                        h_baseline = h_baseline + h_baseline*d2*d_acute;
                        h_baseline_individual = h_baseline;
                      }
                     }
                     
                     if(I_temp == 1){
                       if(t >= t_clear_temp){
                         I_vec[j] = 0;
                       }
                       if(t < t_clear_temp){
                         I_vec[j] = 1;
                       }
                     }
                     

                    //Rprintf(\"infected is %f\\n\",I_vec[j]);
                     t_clear_vec[j] = t_clear_temp;
                     h_vec[j] = h_temp;
                     h_last_infection_vec[j] = h_last_infection_temp;
                     t_last_infection_vec[j] = t_last_infection_temp;
                     h_baseline_individual_vec[j] = h_baseline_individual;
                     q1_vec[j] = q1;
                     q2_vec[j] = q2;
                     q_vec[j] = q;
                     d_acute_vec[j] = d_acute;
                    // age this individual 
                    age_vec[j] = age_vec[j] + dt/365;
                   } // end HH member loop
                     
                  
                ")

## Measurement model ## --------------------------------------------------------------------------------
dmeasure <- Csnippet("
                     double sig = exp(log_sig);
                     double sig_2 = exp(log_sig_2);
                     double h_observed;
                     double h_latent;
                     double this_lik;
                     
                     // set pointers for arrays
                     double *I_vec = &I_m1;
                     double *h_vec = &h_m1;
                     double *t_clear_vec = &t_clear_m1;
                     double *h_obs_vec = &h_obs_m1;
            
                     int N_MEMBERS = n_members;
                     int j;
                     lik = 0;
                     for( j = 0; j < N_MEMBERS; j ++){

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
                          this_lik = log(pnorm(2,h_latent,sig,1,0));
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
                     double sig2 = exp(log_sig_2);
                     // set pointers for arrays
                     double *h_vec = &h_m1;
                     double *h_obs_vec = &h_obs_m1;

                     int N_MEMBERS = n_members;
                     int j;
                     double h_latent;
                     
                     for( j = 0; j < N_MEMBERS; j ++){
                       double rand;
                        
                        if(log_transform_titers == 1){
                          h_latent = log2(h_vec[j]/10) + 2;
                        }
                        if(log_transform_titers == 0){
                          h_latent = h_vec[j];
                        }
    
                        if(h_latent < 2){
                          rand = rnorm(h_latent,sig2);
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
                          h_obs_vec[j] = rand;
                     }
                     
                     ")
