# Lower bound

## Result

There exists a nonincreasing function $\alpha: [0, 1] \to [0, 1]$ such that $\alpha(1) = 0$ and for all $x \in [0, 1]$
$$
\int\limits_{ 0 }^{ x } \frac{ 1 - \alpha( y )}{ 1 - \alpha( x ) } dy + \int\limits_{ x }^{ 1 } \exp\left({ \int\limits_{ 0 }^{ y } \ln \alpha(w) dw}\right) dy 
\quad \in \quad [0.6653, 0.6720] \, .
$$

## Algorithm

1. Initialize $\alpha_0(\cdot)$.

2. Compute $u_0(\cdot)$, by $	u_0(x) := \int\limits_{0}^{ x } 1 - \alpha_0(x) dx$ 

3. Compute $K(\cdot, u_0(\cdot))$, by $K(x, u_0) := 1 - \exp \left( \int\limits_{0}^{x} \ln( 1 - u_0'(t) ) dt \right)$.

4. For $n = 0 ..$

   1. Define $u_{n + 1}$ as the solution of 
      $$
      \left\{ \begin{array}{l r}	(u'(x))^2 K(x, u_n)		- u''(x) u(x) = 0			& ; x \in (0, 1)\\ 	u'(1) = 1 \\	u(0) = 0 \,.	\end{array} \right.
      $$

   2. Compute $\alpha_{n + 1}(\cdot)$, by $\alpha_{n + 1}(x) = 1 - u_{n + 1}'(x)$.

   3. If $\alpha_{n + 1}(\cdot)$ is "good enough", stop.