
%   Decide resolution
N = 450 ;
grid = 0:(1/(N-1)):1 ; 

%   Propouse a function
% State the function
alpha =	[0.5213    0.5144    0.5041    0.4920    0.4789    0.4650    0.4502    0.4342    0.4186    0.4013    0.3825    0.3620    0.3410    0.3189    0.2954    0.2704    0.2440    0.2162    0.1858    0.1538    0.1193    0.0865    0.0492    0.0000] ;
% Detail function for grid used 
m = length(alpha) - 1 ; 
gridProposal = 0:1/m:1 ;
alpha = interp1( gridProposal, alpha, grid, 'spline') ;
% Define proposal of solution of ode
dx = [ grid(2:end) - grid(1:end-1) , grid(end) - grid(end-1) ] ;
un = cumsum( 1 - alpha ) .* dx ;

%   Compute initial K-factor 
[ K ] = Kfactor( alpha , N );

%   Compute time coefficient
[ factor ] = computeCoefficient( grid, un ) ;

% Plot initial:
%   - Porposal of quantiles
%   - Solution of ode
%   - Prophet inequality factor 
%   - Time coefficient for ode

figure ;

subplot(2, 2, 1), plot( grid, alpha, 'b') ;
subplot(2, 2, 1), ylim( [0, 1] ) ;
subplot(2, 2, 1), xlim( [0, 1] ) ;
subplot(2, 2, 1), legend( '\alpha, quantiles', 'Location', 'Best' ) ;
subplot(2, 2, 1), title( [ 'Quantile proposal' ] ) ;
subplot(2, 2, 2), plot( grid, un, 'r') ;
subplot(2, 2, 2), title( ['Initial -ode solution- proposal, resolution: ', num2str(1/N) ] ) ;
subplot(2, 2, 2), ylim( [0, 1] ) ;
subplot(2, 2, 2), xlim( [0, 1] ) ;
subplot(2, 2, 2), legend( 'u_n', 'Location', 'Best' ) ;
subplot(2, 2, 3), plot( grid, K, '*k') ;
subplot(2, 2, 3), xlim( [0, 1] ) ;
subplot(2, 2, 3), legend( 'Prophet factor', 'Location', 'Best' ) ;
subplot(2, 2, 3), title( [ 'Prophet inquality factor : ', num2str(min(K)) ] ) ;
subplot(2, 2, 4), plot( grid, factor, '-g') ;
subplot(2, 2, 4), xlim( [0, 1] ) ;
subplot(2, 2, 4), ylim( [0, 1] ) ;
subplot(2, 2, 4), legend( 'Time coefficient', 'Location', 'Best' ) ;
subplot(2, 2, 4), title( 'Time coefficient for ode' ) ;

%%      Start iterating

%   Define parameters
%N is already define, is the initial resolution.
NumIter = 10 ;
NumInterIter = 1; %How many iterations before increasing resolution N?
NumVisual = 10 ; % How many iterations before ploting a partial result?

for i = 1:NumIter
    for j = 1:NumInterIter
        %   Solve ode
        un = myOdeSolver( factor, N ) ;

        display( 'terminé exitoso' )
        pause(1) ;
        
        %   Compute new time coefficient
        [ factor ] = computeCoefficient( grid, un ) ;
    end
    
    %   Maybe plot?
    if mod(i,NumVisual) == 0
        
        % Derivative on the grid
        alpha = zeros(1, N) ;
        alpha(1) =  1  -  (  (un(2) - un(1))*N  ) ;
        alpha(2:end-1) =  1  -   (  (un(3:end) - un(1:end-2))*N/2  ) ;
        alpha(end) =  1  -  (  (un(end) - un(end-1))*N  ) ;

        %   Compute initial K-factor 
        [ K ] = Kfactor( alpha , N );
        
        %   Compute new time coefficient
        [ factor ] = computeCoefficient( grid, un ) ;

        figure ;
        subplot(2, 2, 1), plot( grid, alpha, 'b') ;
        subplot(2, 2, 1), ylim( [0, 1] ) ;
        subplot(2, 2, 1), xlim( [0, 1] ) ;
        subplot(2, 2, 1), legend( '\alpha, quantiles', 'Location', 'Best' ) ;
        subplot(2, 2, 1), title( [ 'Quantile proposal' ] ) ;
        subplot(2, 2, 2), plot( grid, un, 'r') ;
        subplot(2, 2, 2), title( ['Initial -ode solution- proposal, resolution: ', num2str(1/N) ] ) ;
        subplot(2, 2, 2), ylim( [0, 1] ) ;
        subplot(2, 2, 2), xlim( [0, 1] ) ;
        subplot(2, 2, 2), legend( 'u_n', 'Location', 'Best' ) ;
        subplot(2, 2, 3), plot( grid, K, '*k') ;
        subplot(2, 2, 3), xlim( [0, 1] ) ;
        subplot(2, 2, 3), legend( 'Prophet factor', 'Location', 'Best' ) ;
        subplot(2, 2, 3), title( [ 'Prophet inquality factor : ', num2str(min(K)) ] ) ;
        subplot(2, 2, 4), plot( grid, factor, '-g') ;
        subplot(2, 2, 4), xlim( [0, 1] ) ;
        subplot(2, 2, 4), ylim( [0, 1] ) ;
        subplot(2, 2, 4), legend( 'Time coefficient', 'Location', 'Best' ) ;
        subplot(2, 2, 4), title( 'Time coefficient for ode' ) ;
    end

    %   Increase resolution
    %N = N + 100 ;
    %grid = 0:(1/(N-1)):1 ;
    
end

