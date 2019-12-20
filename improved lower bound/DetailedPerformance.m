function [ K ] = DetailedPerformance( x )
%DetailedPerformance is the performance garantee for pice-wise blind quantile
%   strategies, doing the detailed analysis.
    
    % Initialize variables
    m = length(x) ;
    alpha = x ;

    %   Start computations
    
    % Initial building blocks
    Sumalpha = cumsum(1 - alpha) ; 
    Prodalpha = cumprod(alpha) ;
    
    % Detailed factor
    fac = zeros(1, m) ;
    for k = 1:m
        if k <= m/2
            fac(k) = 1 / ( 1 - k/m * ( 1 - alpha(1) ) ) ;
        else
            fac(k) = 2 / ( 1 + alpha(1) );
        end
    end

    % Computing garantee
    K = zeros(m+1, 1) ;
    
    % j = 1
    K(1) = sum( [1 Prodalpha(1:m-1)].^(1/m) .* (1 - alpha.^(1/m)) ./ -log(alpha) ) ;

    % j = m+1
    K(m+1) = sum( 1 - alpha ) / m ;

    % j in {2, ..., m}
    for j = 2:m
        aux = zeros(m-j+1, 1) ;
        for k = j:m
            aux(k-j+1) = Prodalpha(k-1)^(1/m) * fac(k-1) * (1 - alpha(k)^(1/m))/(-log(alpha(k))) ;
        end
        K(j) = Sumalpha(j-1) / (m*(1-alpha(j))) + sum( aux ) ;
    end
end
