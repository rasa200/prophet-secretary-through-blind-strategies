# Upper bound

## Definitions

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

If the solution $\beta_{ K, \overline{t} }$ does not exist, we say that $\alpha_{ K, \overline{t} }$ does not exist either.

## Result

$$
\sup\limits_{ \substack{ K \in [0, 3] \\ \overline{t} \in [0, 1/3] } } \min \left\{
			1 - \int\limits_{ 0 }^{ 1 } \alpha_{ K, \overline{t} }( s ) ds , 
			\int\limits_{ 0 }^{ 1 } e^{ \int\limits_{ 0 }^{ s } \ln \alpha_{ K, \overline{t} }( w ) dw } ds
			\right\} 
		\leq 0.675
$$

## Algorithm

1. Define function
   $$
   f(K, \overline{t}) = \min \left\{ 1 - \int\limits_{ 0 }^{ 1 } \alpha_{ K, \overline{t} }( s ) ds ,             \int\limits_{ 0 }^{ 1 } e^{ \int\limits_{ 0 }^{ s } \ln \alpha_{ K, \overline{t} }( w ) dw } ds \right\} \,.
   $$

2. Maximize $f(K, \overline{t})$ over $[0, 3] \times [0, 1/3]$ via brute force inspection.

