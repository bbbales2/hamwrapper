data {
  int N;
}

parameters {
  real z[N];
}

model {
  z ~ normal(0, 1);
}
