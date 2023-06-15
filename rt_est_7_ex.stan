//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//



//Example model
data {
  int T; //number of event time
  int S; //number of serial intervals
  int I[T];//observed counts at event time T
  vector[S] w;//serial interval for linked local case
  int tau;//time window
}

transformed data {
  real X[T];
  vector[S] Y;
  for (t in 1:T) {
    for (s in 1:S) {
      if ((t-s)>0)
        Y[s] = I[t-s] * w[s];    
      else 
        Y[s] = 0;
    }
    X[t] = sum(Y);
  }
}


parameters {
  real<lower=0> rt[T];
}

model {
  for (t in 1:T){
    rt[t] ~ gamma(1,0.2); //prior
  }

  for (t in 7:T){
    for (s in (t-tau+1):t){
      if (X[s] > 0)
        I[s] ~ poisson(rt[t]*X[s]);
    }
  }
}



