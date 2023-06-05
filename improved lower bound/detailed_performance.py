# DetailedPerformance is the performance garantee for pice-wise blind quantile strategies, doing the detailed analysis.
import numpy as np

def detailed_performance(alpha):

    # Sanity checks
    for a in alpha:
        if a <= 0 or a >= 1:
            return 0

    # Initialize variables
    m = len(alpha)
    
    #   Start computations
    
    # Detailed factor 
    # g_{m,p}(k), where p = alpha_1
    g = np.zeros(m)
    p = alpha[0]
    for k in range(1, m+1):
        if k <= m/2:
            g[k-1] = 1 / ( 1 - k/m * ( 1 - p ) )
        else:
            g[k-1] = 2 / ( 1 + p )

    print(p, g)

    # Computing garantee
    # f_j( \alpha_1, ..., \alpha_m )
    
    # Initialize
    f = np.ones(m+1)
    
    # j = 1
    f[0] = np.sum( [ np.prod(alpha[0:k-1])**(1/m) * (1 - alpha[k-1]**(1/m)) / -np.log(alpha[k-1])  for k in range(1, m + 1)])

    # j = m+1
    f[m] = np.sum( 1 - alpha ) / m ;

    # j in {2, ..., m}
    # Initial building blocks
    Sumalpha = np.cumsum(1 - alpha) 
    Prodalpha = np.cumprod(alpha)
    # Define for each j
    for j in range(2, m+1):
        aux = [Prodalpha[k-2]**(1/m) * g[k-2] * (1 - alpha[k-1]**(1/m)) / (-np.log(alpha[k-1]))  for k in range(j, m+1)]
        f[j-1] = Sumalpha[j-2] / ( m * (1 - alpha[j-1]) ) + np.sum( aux ) 

    return f
