data{
    int family[87];
    vector[87] value;
    int site[87];
    int habitat_dummy[87];
    int species[87];
    int benthic_pelagic_dummy[87];
    vector[87] tm;
    vector[87] K;
    vector[87] Lmax;
    vector[87] trophic_level;
}
parameters{
    vector[9] B8_sp;
    vector[9] a_sp;
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
    vector[87] mu;
    Rho ~ lkj_corr( 3 );
    scale ~ exponential( 1 );
    sigma_site ~ exponential( 1 );
    sigma_species ~ exponential( 1 );
    X1 ~ normal( 0 , 1 );
    B8 ~ normal( 0 , 1 );
    B7 ~ normal( 0.0752 , 0.119 );
    B6 ~ normal( 0.121 , 0.144 );
    B5 ~ normal( -0.0441 , 0.158 );
    B4 ~ normal( -0.245 , 0.126 );
    B3 ~ normal( -0.208 , 0.124 );
    a ~ normal( -0.131 , 0.172 );
    {
    vector[2] YY[9];
    for ( j in 1:9 ) YY[j] = [ a_sp[j] , B8_sp[j] ]';
    YY ~ multi_normal( rep_vector(0,2) , quad_form_diag(Rho , sigma_species) );
    }
    for ( i in 1:87 ) {
        mu[i] = B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + a + a_sp[species[i]] + (B8 + B8_sp[species[i]]) * habitat_dummy[i] + X1[site[i]] * sigma_site;
        mu[i] = exp(mu[i]);
    }
    value ~ gamma( mu/scale , 1/scale );
}
generated quantities{
    vector[87] log_lik;
    vector[87] mu;
    for ( i in 1:87 ) {
        mu[i] = B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + a + a_sp[species[i]] + (B8 + B8_sp[species[i]]) * habitat_dummy[i] + X1[site[i]] * sigma_site;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:87 ) log_lik[i] = gamma_lpdf( value[i] | mu/scale , 1/scale );
}

