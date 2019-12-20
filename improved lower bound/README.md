# Improved lower bound

## Definitions

Let us define for $m \geq 1$ and $p \in (0, 1)$ the function $g_{ m , p } : [m] \to \mathbb{R}_+$ by 
$$
g_{ m , p }( k ) = \left\{ \begin{array}{l l}
	\frac{ 1 }{ 1 - \frac{k}{m} ( 1 - p ) }
    &; k \leq \frac{m}{2} \\
    \frac{ 2 }{ 1 + p }
    &; k > \frac{m}{2} \,.
    \end{array}	 
    \right.
$$
Given $1 \geq \alpha_1 \geq \ldots \geq \alpha_m > 0$, define $\alpha_{\alpha_1, \ldots, \alpha_m}: [0, 1] \to [0, 1]$ by 
$$
\alpha_{\alpha_1, \ldots, \alpha_m}( x )
		= \sum\limits_{ j \in [m] } \alpha_j \mathbb{1}_{ \left[ \frac{j-1}{m}, \frac{j}{m} \right) }( x ) \,.
$$

## Result

There exists $1 \geq \alpha_1 \geq \ldots \geq \alpha_m > 0$ such that, for every instance $F_1, \ldots, F_mN$ and $t > 0$,  
$$
\mathbb{P}( V_{\sigma_T} > t ) 
	\geq \min_{ j \in [m+1] } \{ f_j( \alpha_1, \ldots, \alpha_m ) \} \cdot \mathbb{P}( \max\limits_{ i \in [n] } \{ V_i \} > t) \,,
$$
with
$$
f_j( \alpha_1, \ldots, \alpha_m ) = \left\{ \begin{array}{l l}
​		\sum\limits_{ k = 1 }^{ m } 
​			\left( \prod\limits_{ l \in [k-1] } \alpha_l \right)^{ \frac{1}{m} } 
​			\left( \frac{ 1 - \alpha_k^{ \frac{1}{m} } }{ -\ln \alpha_k } \right)
​			&; j = 1 \\
​		\sum\limits_{ k \in [j-1] } \frac{ 1 - \alpha_k }{ m( 1 - \alpha_j ) } 

  + \sum\limits_{ k = j }^{ m } 
    	\left( \prod\limits_{ l \in [k-1] } \alpha_l \right)^{ \frac{1}{m} } 
    	g_{ m, \alpha_1 }( k - 1 ) 
    	\left( \frac{ 1 - \alpha_k^{ \frac{1}{m} } }{ -\ln \alpha_k } \right)
    		&; j \in \{ 2, \ldots, m \} \\
    	\sum\limits_{ k \in [ m ] } \frac{ 1 - \alpha_k }{ m }
    		&; j = m+1 \,.
    	\end{array}
    	\right.
$$
where $T$ is the stopping time defined by the blind strategy  $\alpha_{\alpha_1, \ldots, \alpha_m}$.

## Algorithm

1. Take $m = 30$.
2. Find $\alpha^* \in argmax\{ \min_{ j \in [m+1] } \{ f_j( \alpha_1, \ldots, \alpha_m ) \} \}$.