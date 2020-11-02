data{
    vector[145] logvalue;
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
    real<lower=0> sigma_fam;
    real<lower=0> sigma_site;
    real<lower=0> sigma;
}
model{
    vector[145] mu;
    sigma ~ exponential( 2 );
    sigma_site ~ exponential( 2 );
    sigma_fam ~ exponential( 2 );
    B7 ~ normal( -0.159 , 0.253 );
    B6 ~ normal( -0.065 , 0.224 );
    B5 ~ normal( -0.857 , 0.48 );
    B4 ~ normal( -1.014 , 0.207 );
    B3 ~ normal( -0.131 , 0.237 );
    a ~ normal( 8.212 , 1.14 );
    X2 ~ normal( 0 , 1 );
    X1 ~ normal( 0 , 1 );
    for ( i in 1:145 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
    }
    logvalue ~ normal( mu , sigma );
}
generated quantities{
    vector[145] log_lik;
    vector[145] mu;
    for ( i in 1:145 ) {
        mu[i] = a + B3 * trophic_level[i] + B4 * Lmax[i] + B5 * K[i] + B6 * tm[i] + B7 * benthic_pelagic_dummy[i] + X1[site[i]] * sigma_site + X2[family[i]] * sigma_fam;
    }
    for ( i in 1:145 ) log_lik[i] = normal_lpdf( logvalue[i] | mu[i] , sigma );
}

