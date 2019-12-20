function [ Factors ] = TraslationFactors( tgrid, gamma, tbargrid )
%

    dt = tgrid(2:end) - tgrid(1:end-1) ;

    lngamma = cumsum( log(gamma).* dt ) ;
    
    Factors = zeros( 1, length(tbargrid) ) ;

    for i = 1:length(tbargrid)
        tbar = tbargrid(i) ;
        Factors(i) = (1 - tbar) * mean( (1 - gamma(1:end-1)) .* (dt) ) ;
        Factors(i) = min( Factors(i), tbar + (1 -tbar) * sum( exp( (1-tbar) .* lngamma ) .* dt ) ) ;
    end
    
    


end

