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

## Cached results

- alpha = [0.62015846, 0.60632896, 0.60132068, 0.59167718, 0.56560115,
       0.55355759, 0.50561937, 0.44764309, 0.43340237, 0.41983029,
       0.40474392, 0.37860404, 0.36393637, 0.32562933, 0.3123375 ,
       0.29555341, 0.28525343, 0.28213181, 0.26256917, 0.24873521,
       0.23227078, 0.21568684, 0.19692467, 0.17770401, 0.15862058,
       0.13520886, 0.11215676, 0.0875786 , 0.06261711, 0.03116846]
- performance = 0.66951433333

## Implementations

Matlab, Python and Rust implementations are provided. 

The Matlab script gives the best answer so far (performance $0.66975$).

Python (scipy) has given the solution
- alpha = [0.67699182, 0.66490484, 0.60158186, 0.55614734, 0.52010464,
   0.49025544, 0.46484012, 0.4424727 , 0.4224104 , 0.40415579,
   0.3874404 , 0.37196344, 0.35793043, 0.34318325, 0.32944535,
   0.31320542, 0.29739753, 0.28191541, 0.26670487, 0.25143146,
   0.23601253, 0.22021746, 0.20378073, 0.18680571, 0.16867315,
   0.14915083, 0.12761786, 0.10292822, 0.07331355, 0.03333235] 
- performance = 0.6684561705990512

Rust (argmin crate, simulated annealing) has given the solution
- alpha =  [0.6488377139147816, 0.6286748255066426, 0.5959259780730781, 0.581146287337283, 0.5127584720607117, 0.49709499001209534, 0.4800210595809111, 0.46880034959807604, 0.4383713763928919, 0.4073292082536808, 0.36363890895975204, 0.35956533568360116, 0.34557391150818206, 0.34189193603737095, 0.33733738779663586, 0.3292091977620634, 0.3264539322982454, 0.2821942725093516, 0.263229188864773, 0.24683708037169247, 0.24287176137061234, 0.21895174321906052, 0.21534660761148774, 0.1877662629954155, 0.1580734012651116, 0.14218799779672164, 0.11287582062790322, 0.10297556903027627, 0.07025581930535822, 0.04094300669695767]
- performance = 0.6684144764570266
