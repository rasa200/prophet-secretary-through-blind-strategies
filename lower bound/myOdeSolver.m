function [ un ] = myOdeSolver( factor, N )
%   Compute ODE by ode15s: M(t, u) u'(t) = f(t, u)
%
%       |   1   0   |   |  u  |'     _  |  u'(t)               |
%       |   0  u(t) |   |  u' |(t)   -  |  (u'(t))^2 * f(t)    |
%
%       u(0) = 0
%       u'(1) = 1
%
%   INPUT:
%       -- factor -- vector, length = n. Corresponds to f(t), evaluated in
%           the grid 0, 1/n, 2/n, ..., 1.
%       -- N -- int, > 0. Desired resolution of the solution, ie: evaluated
%           in the grid 0, 1/N, 2/N, ..., 1.
%
%   OUTPUT:
%       -- un -- vector, length = N. Numerical solution for the ode,
%           evaluated in the grid 0, 1/N, 2/N, ..., 1.
%
%   METHOD:
%       ode15s, Matlab ode-solver. 
%       Shooting method to get the initial/final conditions.

    %   Flag for exceptions
    flag = 0 ;

    %   Time span
    tspan = [ 1 0 ] ; % From t = 1 to t = 0.
    
    %   Transform the time-coefficient in handle function
    n = length(factor) ;
    factor = @(t) interp1(0:1/(n-1):1, factor, [t], 'spline') ;
    
    %   Solve by shooting method
    
    % State eqation
    %   Right-hand side
    fun = @(t,x) [x(2) ; (x(2)).^2 * factor(t) ] ;
    % Left-hand side mass function
    M = @(t, x) [ 1 0 ; 0 x(1) ] ;
    
    % Resolution of shooting corresponds to desired resolution N
    AchievedShooting = 0.6 ;
    for i = 2:floor(log(N))
        % Find the i-th decimal of shooting
        for decimal = 0:9
            
            %   Initial conditions
            shooting = AchievedShooting + decimal*10^(-i) 
            if and( flag == 1 , decimal > 0 )
                shooting = LastShooting - 10^(-i)
            end
            %       Shooting until: u(0) = 0 and u'(0) > 0, 
            %           ie: u(end, 1) = 0 and u(end,2) > 0.
            %   How to do this?
            %   If u'(0) = 0                ->      increase shooting
            %   If u explodes before zero   ->      decrease shooting
            %   If u(0) > 0                 ->      decrease shooting
            %   If u(0) < 0                 ->      increase shooting

            IniCond = [ shooting ; 1 ] ; 
            options = odeset( 'Mass', M, 'MassSingular', 'yes' , 'MStateDependence', 'weak' , 'InitialSlope', [1 factor(end)/shooting]);

            
            
            %   Compute solution
            try
                [ tgrid, u ] = ode15s(fun, tspan, IniCond, options) ;
            end
            
            
            display( ['llegué hasta: ', num2str( tgrid(end) ) ] )
            
            if tgrid(end) > 0
                % Si no encontré en esta vuelta,
                if decimal == 9
                    % Guardo la última solución.
                    if flag == 0
                        LastShooting = shooting + 10^(-i) ;
                        flag = 1 ;
                    end
                end
                % Si ya no puedo agregar más decimales, 
                if and( i == floor(log(N)), decimal == 9)
                    display( 'llegamos al final, recuperemos solución positiva' )
                    shooting = LastShooting ;
                    IniCond = [ shooting ; 1 ] ; 
                    options = odeset( 'Mass', M, 'MassSingular', 'yes' , 'MStateDependence', 'weak' , 'InitialSlope', [1 factor(end)/shooting]);
                    %   Compute solution
                    try
                        [ tgrid, u ] = ode15s(fun, tspan, IniCond, options) ;
                    end
                end
            end
            
            %   Check conditions for shooting
            if tgrid(end) == 0
                
                if u(end, 1) < 0
                    display( 'llegué, pero soy negativo' )
                    if decimal == 9
                        AchievedShooting = shooting 
                        if i == floor(log(N))
                            display( 'llegamos al final, recuperemos solución positiva' )
                            shooting = shooting + 10^(-i)
                            IniCond = [ shooting ; 1 ] ; 
                            options = odeset( 'Mass', M, 'MassSingular', 'yes' , 'MStateDependence', 'weak' , 'InitialSlope', [1 factor(end)/shooting]);
                            %   Compute solution
                            try
                                [ tgrid, u ] = ode15s(fun, tspan, IniCond, options) ;
                            end
                        end
                        break ;
                    end
                end 
                
                if u(end,1) > 0 
                    display( 'llegué y soy positivo' )
                    %Do we have more decimals?
                    if i < floor(log(N))
                        AchievedShooting = shooting - 10^(-i) 
                        LastShooting = shooting ;
                        break ;
                    else
                        display( 'llegamos al final y es solución positiva ' )
                        AchievedShooting = shooting 
                        break ;
                    end
                end
            end
        end
    end

    %   Transform solution to a new proposal for \alpha
    aux = @(x) interp1( tgrid, u(:, 1), [x], 'spline' ) ;

    un = zeros(1,N) ;
    grid = 0:1/(N-1):1 ;
    
    for i = 1:N
        un(i) = aux( grid(i) ) ;
    end

end

