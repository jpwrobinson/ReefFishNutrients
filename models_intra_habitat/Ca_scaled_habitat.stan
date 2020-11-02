data{
    int family[65];
    vector[65] logvalue;
    int site[65];
    int habitat_dummy[65];
    int species[65];
    int benthic_pelagic_dummy[65];
    vector[65] tm;
    vector[65] K;
    vector[65] Lmax;
    vector[65] trophic_level;
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
    vector[7] X1;
    vector<lower=0>[2] sigma_species;
    real<lower=0> sigma_site;
    real<lower=0> sigma;
    corr_matrix[2] Rho;
}
model{
    vector[65] mu;
    Rho ~ lkj_corr( 3 );
    sigma ~ exponential( 2 );
    sigma_site ~ exponential( 1 );
    sigma_species ~ exponential( 1 );
    X1 ~ normal( 0 , 1 );
    B8 ~ normal( 0 , 1 );
    B7 ~ normal( -0.159 , 0.253 );
    B6 ~ normal( -0.065 , 0.224 );
    B5 ~ normal( -0.857 , 0.48 );
    B4 ~ normal( -1.014 , 0.207 );
    B3 ~ normal( -0.131 , 0.237 );
    a ~ normal( 0 , 100 );
    {
    vector[2] YY[9];
    for ( j in 1:9 ) YY[j] = [ a_sp[j] , B8_sp[j] ]';
    YY ~ multi_normal( rep_vector(0,2) , quad_form_diag(Rho , sigma_species) );
    }
    for ( i in 1:65 ) {
        mu[i] = B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + a + a_sp[species[i]] + (B8 + B8_sp[species[i]]) * habitat_dummy[i] + X1[site[i]] * sigma_site;
    }
    logvalue ~ normal( mu , sigma );
}
generated quantities{
    vector[65] log_lik;
    vector[65] mu;
    for ( i in 1:65 ) {
        mu[i] = B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + a + a_sp[species[i]] + (B8 + B8_sp[species[i]]) * habitat_dummy[i] + X1[site[i]] * sigma_site;
    }
    for ( i in 1:65 ) log_lik[i] = normal_lpdf( logvalue[i] | mu[i] , sigma );
}

