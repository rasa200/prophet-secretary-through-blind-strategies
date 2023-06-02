# prophet-secretary-through-blind-strategies

Collection of scripts that solve computational problems stated in the paper [Prophet Secretary Through Blind Strategies](https://doi.org/10.1137/1.9781611975482.118).

There are three results that need the help of computations to be proven and are listed below.

## Lower bound

### Result

There exists a nonincreasing function $\alpha: [0, 1] \to [0, 1]$ such that 
1. $\alpha(1) = 0$, and
2. For all $x \in [0, 1]$, we have that 

$$
\int_{ 0 }^{ x } \frac{ 1 - \alpha(y) }{ 1 - \alpha(x) } dy + \int_{ x }^{ 1 } \exp \left( \int_{ 0 }^{ y } \ln( \alpha(w) ) dw \right) dy \quad \in \quad [0.6653, 0.6720] .
$$


## Improved lower bound

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

## Upper bound

### Definitions

For $K \in \[0, 3\]$ and $\overline{t} \in \[0, 1/3\]$, define 

$$
\alpha_{ K, \overline{t} }( t ) = \begin{cases}
		1
			&; 0 \leq t < \overline{t} \\
		\beta_{ K, \overline{t} }( t )
			&; \overline{t} \leq t \leq 1,
	\end{cases}
$$

where $\beta_{ K, \overline{t} }$ is the solution of the following integro-differential equation

$$
\begin{align}
	\begin{cases}
		\dot{\beta}( t ) = 
			- K \exp\left\[ \int\limits_{ \overline{t} }^{ t } \ln \beta( s ) ds \right\]
    			& t \in ( \overline{t}, 1) \\
	    	\beta( 1 ) = 0
	\end{cases}
\end{align}
$$

### Result

$$
\sup_{ K \in \[0, 3\] \quad \overline{t} \in \[0, 1/3\] } \min \left( 
		1 - \int\limits_{ 0 }^{ 1 } \alpha_{ K, \overline{t} }( s ) ds , 
		\int\limits_{ 0 }^{ 1 } e^{ \int\limits_{ 0 }^{ s } \ln \alpha_{ K, \overline{t} }( w ) dw } ds
	\right)
	\le 0.675
$$
