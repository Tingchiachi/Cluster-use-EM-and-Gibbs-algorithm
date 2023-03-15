# Cluster-use-EM-and-Gibbs-algorithm
EM(Expectation-maximization) and Gibbs algorithm

## Data
[Old Faithful Geyser Data](https://search.r-project.org/CRAN/refmans/tclust/html/geyser2.html)  
A bivariate data set obtained from the Old Faithful Geyser, containing the eruption length and the length of the previous eruption for 271 eruptions of this geyser in minutes.
## EM(Expectation-maximization) Algorithm
We want to maximize a likelihood function.  
**First,repeating E-step.**  
Given initial $\theta^{(0)} = (\alpha_j,\mu_j,\Sigma_j,j=1,2)$,  
updating $\theta^{(t)}$ to hope $Q(\theta|\theta^{(t)})$ converges.  
**Second,repeating M-step to make $Q(\theta|\theta^{(t)})$ converges.**

<div class="center">
  
|$t$ |   1|   2|   3|   4|   5|
|:--:|:--:|:--:|:--:|:--:|:--:|
|$Q(\theta\mid\theta^{(t)})$|-1222.725|-1157.330|-1138.504|-1131.450|-1130.925|

|$t$ |   6|   7|   8|   9|  10|
|:--:|:--:|:--:|:--:|:--:|:--:|
|$Q(\theta\mid\theta^{(t)})$|-1130.946|-1130.955|-1130.958|-1130.959|-1130.959|

</div>

If we set convergence condition as $|\theta^{(t)}-\theta^{(t-1)}|<10^{-3}$.  
We can find when $t=9$, $Q(\theta|\theta^{(t)})$ converges.  
So,the cluster figure is as following.  
![image](https://github.com/Tingchiachi/Cluster-use-EM-and-Gibbs-algorithm/blob/main/em.jpeg)

## Gibbs Algorithm
First,given initial $\theta^{(0)} = (\alpha_j,\mu_j,\tau_j=\sigma_j^{-2},j=1,2)$,  
updating $\theta^{(t)}$ to hope $\theta^{(t)}$ converges.We consider full conditionals:  
$$p(y_i^{(t+1)}|\mu^{(t)},\tau^{(t)},\alpha^{(t)},y^{(t)}_{-i},x)\sim categorical(p_1,p_2),i=1,2,\dots,n$$

$$p(\mu_j^{(t+1)}|\mu^{(t)}_{(-j)},\tau^{(t)},\alpha^{(t)},y^{(t)},x)\sim \mathcal{N}(\bar{\mu_j},\bar{\tau^{-1}_j}),j=1,2$$

$$p(\tau^{(t+1)}_j|\mu^{(t)},\tau^{(t)}_{(-j)},\alpha^{(t)},y^{(t)},x)\sim Gamma(a+\dfrac{n_j}{2},b+\dfrac{1}{2}\sum_{i;y_i=j}(x_i-\mu_j)^2),j=1,2$$

$$p(\alpha^{(t+1)}|\mu^{(t)},\tau^{(t)},\alpha^{(t)}_{(-j)},y^{(t)},x)\sim Dirichlet(a_1+n_1,a_2+n_2)$$

Iterating 1000 times,the cluster figure is as following.
![image](https://github.com/Tingchiachi/Cluster-use-EM-and-Gibbs-algorithm/blob/main/gibbs.jpeg)
