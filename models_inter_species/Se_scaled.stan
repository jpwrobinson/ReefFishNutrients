data{
    vector[176] value;
    int family[176];
    int site[176];
    int benthic_pelagic_dummy[176];
    vector[176] tm;
    vector[176] K;
    vector[176] Lmax;
    vector[176] trophic_level;
}
parameters{
    vector[8] X1;
    vector[9] X2;
    real a;
    real B3;
    real B4;
    real B5;
    real B6;
    real B7;
    real<lower=0> scale;
    real<lower=0> sigma_fam;
    real<lower=0> sigma_site;
}
model{
    vector[176] mu;
    sigma_site ~ exponential( 2 );
    sigma_fam ~ exponential( 2 );
    scale ~ exponential( 2 );
    B7 ~ normal( 0.524 , 0.278 );
    B6 ~ normal( -0.352 , 0.278 );
    B5 ~ normal( -0.434 , 0.404 );
    B4 ~ normal( 0.441 , 0.192 );
    B3 ~ normal( -0.083 , 0.248 );
    a ~ normal( -0.061 , 1.042 );
    X2 ~ normal( 0 , 0.1 );
    X1 ~ normal( 0 , 1 );
    for ( i in 1:176 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    value ~ gamma( mu/scale , 1/scale );
}
generated quantities{
    vector[176] log_lik;
    vector[176] mu;
    for ( i in 1:176 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:176 ) log_lik[i] = gamma_lpdf( value[i] | mu/scale , 1/scale );
}

