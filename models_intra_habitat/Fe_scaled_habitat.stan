data{
    int family[105];
    vector[105] value;
    int site[105];
    int habitat_dummy[105];
    int species[105];
    int benthic_pelagic_dummy[105];
    vector[105] tm;
    vector[105] K;
    vector[105] Lmax;
    vector[105] trophic_level;
}
parameters{
    vector[10] B8_sp;
    vector[10] a_sp;
    real a;
    real B3;
    real B4;
    real B5;
    real B6;
    real B7;
    real B8;
    vector[8] X1;
    vector<lower=0>[2] sigma_species;
    real<lower=0> sigma_site;
    real<lower=0> scale;
    corr_matrix[2] Rho;
}
model{
    vector[105] mu;
    Rho ~ lkj_corr( 3 );
    scale ~ exponential( 1 );
    sigma_site ~ exponential( 1 );
    sigma_species ~ exponential( 1 );
    X1 ~ normal( 0 , 1 );
    B8 ~ normal( 0 , 1 );
    B7 ~ normal( -0.018 , 0.138 );
    B6 ~ normal( 0.076 , 0.114 );
    B5 ~ normal( 0.034 , 0.229 );
    B4 ~ normal( 0.025 , 0.125 );
    B3 ~ normal( -0.101 , 0.125 );
    a ~ normal( 0 , 100 );
    {
    vector[2] YY[10];
    for ( j in 1:10 ) YY[j] = [ a_sp[j] , B8_sp[j] ]';
    YY ~ multi_normal( rep_vector(0,2) , quad_form_diag(Rho , sigma_species) );
    }
    for ( i in 1:105 ) {
        mu[i] = B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + a + a_sp[species[i]] + (B8 + B8_sp[species[i]]) * habitat_dummy[i] + X1[site[i]] * sigma_site;
        mu[i] = exp(mu[i]);
    }
    value ~ gamma( mu/scale , 1/scale );
}
generated quantities{
    vector[105] log_lik;
    vector[105] mu;
    for ( i in 1:105 ) {
        mu[i] = B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + a + a_sp[species[i]] + (B8 + B8_sp[species[i]]) * habitat_dummy[i] + X1[site[i]] * sigma_site;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:105 ) log_lik[i] = gamma_lpdf( value[i] | mu/scale , 1/scale );
}

