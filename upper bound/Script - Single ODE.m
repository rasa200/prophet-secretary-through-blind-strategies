%   Define factors
tbar = 0.2 ;
alpha_0 = 0.0176;
K = 0.1 ;

%   Compute ODE by ode15s: M(t, u) u'(t) = g(t, u)
%
%       |   1     0       |   |  u  |'     _  |  u'(t)           |
%       |   0  exp[u'(t)] |   |  u' |(t)   -  |  -K exp[u(t)]    |
%
%       u( tbar ) = 0
%       u'( tbar ) = ln( alpha_0 )

%   Time span
tspan = [ tbar 1 ] ;

%       State equation
%   Right-hand side
fun = @(t,x) [ x(2) ; -K * exp( x(1) ) ] ;
%   Left-hand side mass function
M = @(t, x) [ 1 0 ; 0 exp( x(2) ) ] ;
%   Initial conditions
IniCond = [ 0 ; log(alpha_0) ] ; 

%   Options for solver
options = odeset( 'Mass', M, 'MassSingular', 'yes' , 'MStateDependence', 'weak' );

%   Compute solution
[ tgrid, u ] = ode15s(fun, tspan, IniCond, options) ;

exp(u( end, 2 ) )
tgrid(end)

%%
%       Ploting solution

figure ;
plot( tgrid, exp( u(:, 2) ) ) ;

%%
%       Checking conditions over K
tgrid = [ 0:1/100:tbar , tgrid(2:end)' ] ;
alpha = exp( u(:, 2) ) ;
alpha = [ ones( 1 , length(0:1/100:tbar)-1 ) alpha' ] ;

intAlpha = sum( (tgrid( [ 2:end, end] ) - tgrid ) .* alpha ) ;
ntbar = length(0:1/100:tgrid(1)) ;

display( [ 'K = ', num2str(K), ', suspously is = ', num2str( alpha( ntbar ) / ( 1 - intAlpha - tbar ) ) ] ) ;
display( [ 'factor = ', num2str( 1 - intAlpha ) ] ) ;


%%      Working a single K

%   Define parameters
tol = 1e-6 ;
K = 1.2 ;

%   Solve by shooting
[ tgrid, alpha, g ] = myOdeSolver( K, tol ) ;

%%

%   Define translations
tbargrid =  

%   Compute translations factors
[ Factors ] = TraslationFactors( tgrid, gamma, tbargrid )




