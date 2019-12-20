function [ K ] = Kfactor( alpha , m )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    m = m-1;

    % Checking out K
    K = zeros(1, m+1) ;
    % Pre-compute partial products
    parProd = zeros(1, m+1) ;
    for j = 1:m+1
        parProd(j) = prod( alpha(1:j-1) ) ^ (1/m) ;
    end
    % Compute other values

    for j = 1:m+1
        aux = ( 1 - (alpha).^(1/m) ).*parProd ./ (-log(alpha)) ;
        K(j) = sum( 1 - alpha(1:j-1) ) / ( m * (1 - alpha(j)) ) + sum( aux(j:m) ) ;
    end

end

