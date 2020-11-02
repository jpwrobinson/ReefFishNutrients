data{
    vector[150] value;
    int family[150];
    int site[150];
    int benthic_pelagic_dummy[150];
    vector[150] tm;
    vector[150] K;
    vector[150] Lmax;
    vector[150] trophic_level;
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
    vector[150] mu;
    sigma_site ~ exponential( 2 );
    sigma_fam ~ exponential( 2 );
    scale ~ exponential( 2 );
    B7 ~ normal( 0.008 , 0.114 );
    B6 ~ normal( 0.123 , 0.095 );
    B5 ~ normal( 0.032 , 0.19 );
    B4 ~ normal( -0.155 , 0.867 );
    B3 ~ normal( -0.194 , 0.091 );
    a ~ normal( -0.399 , 0.418 );
    X2 ~ normal( 0 , 1 );
    X1 ~ normal( 0 , 1 );
    for ( i in 1:150 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    value ~ gamma( mu/scale , 1/scale );
}
generated quantities{
    vector[150] log_lik;
    vector[150] mu;
    for ( i in 1:150 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:150 ) log_lik[i] = gamma_lpdf( value[i] | mu/scale , 1/scale );
}

