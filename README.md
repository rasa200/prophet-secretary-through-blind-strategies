# prophet-secretary-through-blind-strategies

Collection of scripts that solve computational problems stated in the paper "Prophet secretary through blind strategies".



There are three results that need the help of computations to be proven and are listed below.

## Lower bound

### Result

There exists a nonincreasing function $\alpha: [0, 1] \to [0, 1]$ such that $\alpha(1) = 0$ and for all $x \in [0, 1]$
$$
\int\limits_{ 0 }^{ x } \frac{ 1 - \alpha( y )}{ 1 - \alpha( x ) } dy + \int\limits_{ x }^{ 1 } \exp\left({ \int\limits_{ 0 }^{ y } \ln \alpha(w) dw}\right) dy 
\quad \in \quad [0.6653, 0.6720] \, .
$$

### Script



## Improved lower bound

### Definitions

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

### Result

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



## Upper bound for blind strategies

### Definitions

For $K \in [0, 3]$ and $\overline{t} \in [0, 1/3]$, define 
$$
\alpha_{ K, \overline{t} }( t ) = \left\{ \begin{array}{l l}
		1
			&; 0 \leq t < \overline{t} \\
		\beta_{ K, \overline{t} }( t )
			&; \overline{t} \leq t \leq 1,
		\end{array}
		\right.
$$
where $\beta_{ K, \overline{t} }$ is the solution of the following integro-differential equation
$$
\begin{align}
	\left\{ \begin{array}{l l}
		\dot{\beta}( t ) = 

  - K \exp\left[ \int\limits_{ \overline{t} }^{ t } \ln \beta( s ) ds \right]
    		& t \in ( \overline{t}, 1) \\
    	\beta( 1 ) = 0 \,.
    	\end{array}
    	\right.
    \end{align}
$$

### Result

$$
\sup\limits_{ \substack{ K \in [0, 3] \\ \overline{t} \in [0, 1/3] } } \left\{
			1 - \int\limits_{ 0 }^{ 1 } \alpha_{ K, \overline{t} }( s ) ds , 
			\int\limits_{ 0 }^{ 1 } e^{ \int\limits_{ 0 }^{ s } \ln \alpha_{ K, \overline{t} }( w ) dw } ds
			\right\} 
		\leq 0.675
$$

