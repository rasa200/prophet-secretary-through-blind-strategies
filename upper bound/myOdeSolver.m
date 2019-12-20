function [ tgrid, alpha, g ] = myOdeSolver( K, tbar, tol )
% Compute ODE by ode15s: M(t, u) u'(t) = g(t, u)
%
%       |   1     0       |   |  g  |'     _  |  g'(t)           |
%       |   0  exp[u'(t)] |   |  g' |(t)   -  |  -K exp[g(t)]    |
%
%       g( tbar ) = 0
%       g'( 1 ) = -oo
%
% Uses the shooting method in g'(0) = ln( alpha_0 ).
% INPUT:
%   - K - double. Defines ODE. It should be in [1, 3/2].
%   - tol - double. Tolerance for the shooting method. It allows that g'(0)
%       might be tol wrong.
% OUTPUT:
%   - tgrid - list. Times where the solution was computed.
%   - alpha - list. Corresponds to exp[ g'(.) ] at tgrid.
    
    %   Define error control parameters
    RelTol = 1e-10 ; % odeget(odeset, 'RelTol') ; % 1e-8 ;
    AbsTol = 1e-10 ; % odeget(odeset, 'AbsTol') ; % 1e-8 ;
    
    Refine = 1000 ; % How many refinement for each step?
    MaxStep = 0.1 ; % Last resource!!! (prefere Refine) max step distance

    %   Time span
    tspan = [ tbar 1 ] ;

    %       State equation
    %   Right-hand side
    fun = @(t,x) [ x(2) ; -K * exp( x(1) ) ] ;
    %   Left-hand side mass function
    M = @(t, x) [ 1 0 ; 0 exp( x(2) ) ] ;
    %   Options for solver
    options = odeset( 'Mass', M, 'MassSingular', 'yes' , 'MStateDependence', 'weak', 'RelTol', RelTol, 'AbsTol', AbsTol, 'Refine', Refine, 'MaxStep', MaxStep );
    
    
    %   Start shooting method
    % We know that g'(0) = ln( alpha_0 ), 
    % with alpha_0 in [1/2, 1].
    down = 0.002 ;
    up = 1 ;
    
    while up - down > tol
        %       Compute solution at middle point
        alpha_0 = mean([ down, up ]) ;
        
        %   Initial conditions
        IniCond = [ 0 ; log(alpha_0) ] ; 
        %   Compute solution
        [ tgrid, g ] = ode15s(fun, tspan, IniCond, options) ;
        
        %       Clasify middle point
        if tgrid(end) == 1
            up = alpha_0 ;
        else
            down = alpha_0 ;
        end
    end

    alpha = exp( g(:, 2) ) ;
    g = g(:, 1) ;

end

