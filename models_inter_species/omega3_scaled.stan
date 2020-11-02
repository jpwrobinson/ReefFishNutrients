data{
    vector[145] value;
    int family[145];
    int site[145];
    int benthic_pelagic_dummy[145];
    vector[145] tm;
    vector[145] K;
    vector[145] Lmax;
    vector[145] trophic_level;
}
parameters{
    vector[7] X1;
    vector[10] X2;
    real a;
    real B3;
    real B4;
    real B5;
    real B6;
    real B7;
    real scale;
    real<lower=0> sigma_fam;
    real<lower=0> sigma_site;
}
model{
    vector[145] mu;
    sigma_site ~ exponential( 1 );
    sigma_fam ~ exponential( 1 );
    scale ~ cauchy( 0 , 0.5 );
    B7 ~ normal( -0.032 , 0.189 );
    B6 ~ normal( -0.108 , 0.183 );
    B5 ~ normal( -0.14 , 0.436 );
    B4 ~ normal( -0.174 , 0.12 );
    B3 ~ normal( 0.013 , 0.165 );
    a ~ normal( 2.354 , 0.84 );
    X2 ~ normal( 0 , 1 );
    X1 ~ normal( 0 , 1 );
    for ( i in 1:145 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    value ~ gamma( mu/scale , 1/scale );
}
generated quantities{
    vector[145] log_lik;
    vector[145] mu;
    for ( i in 1:145 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:145 ) log_lik[i] = gamma_lpdf( value[i] | mu/scale , 1/scale );
}

