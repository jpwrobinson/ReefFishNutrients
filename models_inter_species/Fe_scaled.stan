data{
    vector[179] value;
    int family[179];
    int site[179];
    int benthic_pelagic_dummy[179];
    vector[179] tm;
    vector[179] K;
    vector[179] Lmax;
    vector[179] trophic_level;
}
parameters{
    vector[8] X1;
    vector[10] X2;
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
    vector[179] mu;
    sigma_site ~ exponential( 2 );
    sigma_fam ~ exponential( 2 );
    scale ~ exponential( 2 );
    B7 ~ normal( -0.018 , 0.138 );
    B6 ~ normal( 0.076 , 0.114 );
    B5 ~ normal( 0.034 , 0.229 );
    B4 ~ normal( 0.025 , 0.125 );
    B3 ~ normal( -0.101 , 0.125 );
    a ~ normal( -1.875 , 0.645 );
    X2 ~ normal( 0 , 1 );
    X1 ~ normal( 0 , 1 );
    for ( i in 1:179 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    value ~ gamma( mu/scale , 1/scale );
}
generated quantities{
    vector[179] log_lik;
    vector[179] mu;
    for ( i in 1:179 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:179 ) log_lik[i] = gamma_lpdf( value[i] | mu/scale , 1/scale );
}

