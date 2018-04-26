Algorithm
-------
NetII-SLDA: Network incorporated integrative sparse linear discriminant analysis


Maintainer
-------
Xiaoyan Wang   <xywang@hnu.edu.cn>


Publication
-------
Network incorporated integrative sparse linear discriminant analysis (2018)


Usage
-------
To apply the algorithm, first use *afunN1* to construct a network structure among predictors, then iteratively update parameters using *update.beta* with given tuning parameters. Sample data can be generated using *gen.data* for experiments.


 
 
1.gen.data

- Usage
 - gen.data(n,beta,sigma,p): generate a single dataset
      
- Arguments
 - n: sample size
 - beta: length-p vector. Represent the true coefficients of subgroup Y=2
 - sigma: a matrix of dimension p. Denote the covariance matrix among predictors X
 - p: the dimension of X
      
- Values
 - This function returns a list containing six objects named y, y1, y2, x1, x2, and x.sub. Here y is a length-n response observation vector and x.sub is a  design matrix.  “1” and “2” denote the class labels, so y1 and x1 are the observations with class label y=1.

2.	afunN1. This function is defined based on the N.1 adjacency measure.

- Usage 
  - afunN1(alpha,x.org): construct the N1 type adjacency matrix for observation x.org
      
- Arguments
  - alpha: parameter of the N.1 adjacency measure
  - x.org: a matrix of dimension nobs  nvars
      
- Values
  - It returns an adjacency matrix of dimension nvars  nvars

3.	update.beta

- Usage
  - update.beta (beta0,lambda1,lambda2,gamma): update the parameters  given tunings and initial values
      
- Arguments
  - beta0: a vector of the same length with the parameters vector. It is used as an initial value for updating parameters. Suggested choices include OLS, Lasso-estimation, and Ridge-estimation.
  - lambda1: tuning parameter for the 1-norm Group MCP
  - lambda2: tuning parameter for the Laplacian penalty
  - gamma: the tuning of MCP
      
- Values
  - It returns updated parameters of the same dimension as beta0.
