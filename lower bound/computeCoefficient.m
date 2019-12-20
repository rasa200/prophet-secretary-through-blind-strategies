function [ factor ] = computeCoefficient( grid, f )
% computeFactor computes the following expresion
%   1 - exp [ \int_{0}^{t} ln( 1 - f'(x) ) dx ].
%   INPUT:
%       -- grid -- vector. Domain and points where f is known.
%       -- f -- vector. Value of f over grid.

    %   Initialie parameters
    N = length(grid) ;
    f1 = zeros(1, N) ;
    dx = [ grid(2:end) - grid(1:end-1) , grid(end) - grid(end-1) ] ;
    
    %   Compute derivative.
    % First point
    f1(1) = ( f(2) - f(1) ) ;
    % Middle points
    f1(2:end-1) = ( f(3:end) - f(1:end-2) ) ./ 2 ;
    % End point
    f1(end) = ( f(end) - f(end-1) ) ;
    % Re-scale
    f1 = f1 ./ dx ;

    %   Compute intengrand
    integrand = log( 1 - f1 ) ;
    
    %   Compute integral
    I = [ cumsum( integrand ) .* dx ] ;
    
    %   Compute factor
    factor = 1 - exp(I) ;
    
    %   Correct last term
    factor(end) = interp1( grid(1:end-1), factor(1:end-1), [1], 'spline' ) ;
end



