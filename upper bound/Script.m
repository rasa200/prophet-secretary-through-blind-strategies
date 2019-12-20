%   Define parameters
Ks = 0.1:0.5:6 ;
Ks = 0.1:0.2:2 ;
Ks = 1:0.1:1.6 ;
Ks = [ 1.05:0.02:1.15  ,  1.17:0.02:1.25 ] ;
Ks = [ 1:0.2:1.5 1.1:0.05:1.3 1.15:0.01:1.25 1.20:0.004:1.24 ];
Ks = 1.1:0.02:1.3 ;
Ks = 1.15:0.01:1.25 ;
Ks = 1.20:0.004:1.24 ;
Ks = 1.21:0.002:1.23 ;
Ks = 1.215:0.001:1.225 ;
Ks = 1.215:0.001:1.218 ;

tbars = 0:0.05:0.5 ;
tbars = 0:0.01:0.2 ;
tbars = [ [ 0 0.01 ] 0.1:0.01:0.14 ] 
tbars = [ [ 0 ] 0.11:0.005:0.15 ] ;
%%
%   Define tolerance values
tol = 1e-6 ; % Shooting tolerance
tolTime = 1e-4 ; % For time grid tolerance, if not small enough...
tolInt = 1e-5 ; % For integrating numerical solution
% For more, see myOdeSolver!!

%   Storing variables
Disparity = zeros( length(Ks), length(tbars)) ;
factorSolution = zeros( length(Ks), length(tbars), 2) ;

for i = 1:length(Ks)
    %   Define parameters
    K = Ks(i) ;
    for j = 1:length(tbars)
        tbar = tbars(j) ;
        %   Compute solution
        try 
            [ tgrid, beta, g ] = myOdeSolver( K, tbar, tol ) ;
            display([ 'Resolví para K = ', num2str(K), ' y tbar = ', num2str(tbar) ]) ;    

            %   Compute interest values
            intBeta = sum( beta .* (tgrid([2:end end]) - tgrid ) ) ;
            PowerInt = sum( exp(g) .* ( tgrid([2:end end]) - tgrid ) ) ;
            factorSolution( i, j, 1 ) = 1 - tbar - intBeta ;
            factorSolution( i, j, 2 ) = tbar + beta(1)/K ;
            Disparity(i, j) = abs( factorSolution( i, j, 1 ) - factorSolution( i, j, 2 ) ) ;
        catch
            display([ 'No existe para K = ', num2str(K), ' y tbar = ', num2str(tbar) ]) ;    
            factorSolution( i, j, 1 ) = NaN ;
            factorSolution( i, j, 2 ) = NaN ;
            Disparity(i, j) = NaN ;
        end
                
    end
end


display('Terminé!!') ;

%%
%       Ploting results
[ X, Y ] = meshgrid( Ks, tbars ) ;

figure ;

subplot(1, 3, 1), surf( X, Y, Disparity') ;
subplot(1, 3, 1), legend( 'Disparity of \alpha_{ K, tbar}', 'Location', 'best' ) ;
subplot(1, 3, 1), xlabel( 'K' ) ;
subplot(1, 3, 1), ylabel( 'tbar' ) ;
subplot(1, 3, 1), title( 'Diference between instances' ) ;

subplot(1, 3, 2), surf( X, Y, factorSolution( :, :, 1 )' ) ;
subplot(1, 3, 2), legend( '\int 1 - \alpha', 'Location', 'best' ) ;
subplot(1, 3, 2), xlabel( 'K' ) ;
subplot(1, 3, 2), ylabel( 'tbar' ) ;
subplot(1, 3, 2), title( 'One variable' ) ;

subplot(1, 3, 3), surf( X, Y, factorSolution( :, :, 2 )' ) ;
subplot(1, 3, 3), legend( '\int exp[\int ln \alpha]', 'Location', 'best' ) ;
subplot(1, 3, 3), xlabel( 'K' ) ;
subplot(1, 3, 3), ylabel( 'tbar' ) ;
subplot(1, 3, 3), title( 'i.i.d. variables' ) ;


figure;

surf( X, Y, min( factorSolution( :, :, 1 )' , factorSolution( :, :, 2 )' ) );
xlabel( 'K' ) ;
ylabel( 'tbar' ) ;
title( 'Overall factor' ) ;