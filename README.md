### R/corihw: Independent Hypothesis Weighting for Correlated Tests


[Yair Goldberg](https://yairgo.net.technion.ac.il/)

[R/corihw](https://github.com/yairgoldy/corihw) is an [R](https:/www.r-project.org) package. This package generalizes [independent hypothesis weighting (IHW)](https://bioconductor.org/packages/release/bioc/html/IHW.html)  multiple testing procedure to correlated tests. It enables using of both correlation between the tests and covarites that are informative under the alternative to increases power over the Bonfferoni correction. 

#### Installation

You can install it from its [GitHub repository](https://github.com/yairgoldy/corihw). You first need to install the [devtools](https://github.com/hadley/devtools) package and [IHW](http://bioconductor.org/packages/release/bioc/html/IHW.html) package.

```r
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IHW")
```

Then install R/corihw using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package.

```r
library(devtools)
install_github("yairgoldy/corihw")
```

#### Example use

Create a 20,000 two-sample t-tests with correlations . Let $\mu=1$ be the mean under the alternative and $\mu=0$ the mean under the null. The proportion of alternatives is 0.05. Let $H_{m\times 1}\in\{0,1\}$ be a vector such that $H_j\sim \text{Bernoulli}(0.05)$ are independent. The first sample of
$15$ observations is drawn from a multivariate random vector $N(H\mu,R)$, and the second sample of $15$ observations is drawn from  $N(0_{m\times 1},R)$. The correlation matrix $R$ is block diagonal with a given block size of 100 and with off-diagonal entries of $\rho=0.9$. 


```{r}

library(corihw)
dat <- toy_example(m=20000,n=15,mu=1,prob=0.05, block_size=100,rho=0.9)

```

Compare IHW, CorIHW, and CorIHW based on M-effective.
```{r}
sol <- corihw (pvalues = dat$P,covariates = dat$X,alpha = 0.1,lambda=5, Sigma=dat$Sigma,
               methods= c("IHW","CorIHW","M-effective"))


rejections <- data.frame(IHW  = sum(sol$rjs_IHW),
                        CorIHW = sum(sol$rjs_Cor),
                        meffective = sum(sol$rjs_meffective))
rejections


```



