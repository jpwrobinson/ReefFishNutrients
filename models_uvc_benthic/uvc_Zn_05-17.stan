data{
    vector[72] value_corr;
    int Year[72];
    int Location[72];
    vector[72] PC1;
}
parameters{
    vector[12] X1;
    vector[6] X2;
    real a;
    real b;
    real<lower=0> sigma_yr;
    real<lower=0> sigma_site;
    real<lower=0> sigma;
}
model{
    vector[72] mu;
    sigma ~ exponential( 2 );
    sigma_site ~ exponential( 2 );
    sigma_yr ~ exponential( 2 );
    b ~ normal( 0 , 10 );
    a ~ normal( 0 , 100 );
    X2 ~ normal( 0 , 1 );
    X1 ~ normal( 0 , 1 );
    for ( i in 1:72 ) {
        mu[i] = a + b * PC1[i] + X1[Location[i]] * sigma_site + X2[Year[i]] * sigma_yr;
    }
    value_corr ~ normal( mu , sigma );
}

