# Cluster-use-EM-and-Gibbs-algorithm
EM(Expectation-maximization) and Gibbs algorithm

## Data

## EM(Expectation-maximization) algorithm
We want to maximize a likelihood function.  
**First,repeating E-step.**  
Given initial $\theta^{(0)} = (\alpha_j,\mu_j,\Sigma_j,j=1,2)$,  
updating $\theta^{(t)}$ to hope $Q(\theta|\theta^{(t)})$ converges.  
**Second,repeating M-step to make $Q(\theta|\theta^{(t)})$ converges.**

\begin{table}[H]
\centering
\caption{Iteration table}
\begin{tabular}{l|rrrrr}
$t$     & 1     & 2     & 3     & 4     & 5 \\
\midrule
$Q(\theta|\theta^{(t)})$   & -1222.725 & -1157.330 & -1138.504 & -1131.450 & -1130.925 \\
\multicolumn{1}{r}{} &       &       &       &       &  \\
$t$     & 6     & 7     & 8     & 9     & 10 \\
\midrule
$Q(\theta|\theta^{(t)})$   & -1130.946 & -1130.955 & -1130.958 & -1130.959 & -1130.959 \\
\end{tabular}%
\label{tab:addlabel}%
\end{table}%

We can find when $t=9$, $Q(\theta|\theta^{(t)})$ converges.  
So,the cluster figure is as following.
![image]()
