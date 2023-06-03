%       Solution for finite m

m = 30;

% First guess
x0 = 0.55:-(0.5/(m-1)):0.05;
%x0 = 0.5.*ones(1, m); % Alternative first guess,
%x0 = x.*0.9 ;

% Linear inequality Ax <= b
A = diag( -ones(m,1) ) + diag( ones(m-1,1), 1 ); % Decreasing alpha
b = zeros(m, 1);

% Linear equality Aeq x = beq
Aeq = [];
beq = [];

% Lower bounds and upper bounds in individual variables
lb = eps*ones(m, 1); % We need to take the log.
ub = ones(m, 1);

% Set options
options = optimoptions('fmincon');
options.Display = 'iter';
options.MaxIter = 1000;
options.TolFun = 1.0e-15;
options.TolCon = 1.0e-15;
options.TolX = 1.0e-50;
options.MaxFunEvals = 100000;

%   Compute solution
x = fmincon(@(x) -min(DetailedPerformance(x)) ,x0,A,b,Aeq,beq,lb,ub,[], options) ;


[ K ] = DetailedPerformance(x) ;
alpha = [x, 0] ;

% Ploting
figure ;
plot(0:1/m:1, alpha) ;
hold on;
plot(0:1/m:1, K) ;
ylim( [0, 1] ) ;
legend( '\alpha', 'K' ) ;
title( [ 'For m = ', num2str(m), ' is K = ', num2str(min(K)) ] ) ;

