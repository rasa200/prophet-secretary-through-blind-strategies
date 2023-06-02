# Improved lower bound

### Definitions

Let us define for $m \geq 1$ and $p \in (0, 1)$ the function $g_{ m , p } : [m] \to \mathbb{R}_+$ by 

$$
g_{m, p}(k) = \begin{cases}
		\frac{ 1 }{ 1 - \frac{k}{m} ( 1 - p ) } &; k \leq \frac{m}{2} \\
		\frac{ 2 }{ 1 + p } &; k > \frac{m}{2}
	\end{cases}
$$

### Result

There exists $1 \ge \alpha_1 \ge \ldots \ge \alpha_m \ge 0$ such that

$$
\min_{j \in [m+1]} f_j( \alpha_1, \ldots, \alpha_m ) \ge 0.66975 \,,
$$

where, for every $j \in [m+1]$,
$$
f_j( \alpha_1, \ldots, \alpha_m ) := \begin{cases}
		\sum\limits_{ k = 1 }^{ m } 
			\left( \prod\limits_{ l \in [k-1] } \alpha_l \right)^{ \frac{1}{m} } 
			\left( \frac{ 1 - \alpha_k^{ \frac{1}{m} } }{ -\ln \alpha_k } \right)
			&; j = 1 \\
		\sum\limits_{ k \in [j-1] } \frac{ 1 - \alpha_k }{ m( 1 - \alpha_j ) } 
			+ \sum\limits_{ k = j }^{ m } 
			\left( \prod\limits_{ l \in [k-1] } \alpha_l \right)^{ \frac{1}{m} } 
			g_{ m, \alpha_1 }( k - 1 ) 
			\left( \frac{ 1 - \alpha_k^{ \frac{1}{m} } }{ -\ln \alpha_k } \right)
				&; j \in \{ 2, \ldots, m \} \\
			\sum\limits_{ k \in [ m ] } \frac{ 1 - \alpha_k }{ m }
				&; j = m+1 \,.
    	\end{cases}
$$

## Algorithm

1. Take $m = 30$.
2. Find $\alpha^* \in argmax\{ \min_{ j \in [m+1] } \{ f_j( \alpha_1, \ldots, \alpha_m ) \} \}$.
