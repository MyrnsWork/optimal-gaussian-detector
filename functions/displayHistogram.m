function fig = displayHistogram( imagetteAmplitude_lin,...
                                 imagettePuissance_lin,...
                                 Pbth_lin,...
                                 iFigure                  )
                             
                             
                             
   fig = figure( iFigure );
   set(fig, 'Units', 'Normalized', 'Position', [0 0 1 1]);

   %----------------------------------------------------------------------
   subplot 121
   N    = 1000;
   Pmax = 1 - 1/length( imagetteAmplitude_lin(:) );
   xMax = raylinv( Pmax, sqrt(Pbth_lin/2) );
   x    = (0 : N-1) * xMax/N;                                  
   fx   = raylpdf( x, sqrt(Pbth_lin/2) );                                   
                                     
   histogram( imagetteAmplitude_lin(:), 'Normalization', 'pdf'), hold on;
   plot( x, fx, 'LineWidth', 1.5), hold off;
   xlabel( "Amplitude" );
   ylabel( "Densité de probabilité" );
   legend("Histogramme", "Densité théorique")
   title( "AMPLITUDE" ); 
   subtitle( "Densité de probabilité avec P_{bth}=" + string( round( 10*log10(Pbth_lin) ) ) + " dB" );  
   grid on;

   %----------------------------------------------------------------------
   subplot 122
   N    = 1000; 
   Pmax = 1 - 1/length( imagettePuissance_lin(:) );
   xMax =  -Pbth_lin * log( 1 - Pmax ) ;
   x    = (0 : N-1) * xMax/N;                                  
   fx   = 1/Pbth_lin * exp( -x/Pbth_lin );                                   
                                     
   histogram( imagettePuissance_lin(:), 'Normalization', 'pdf'), hold on;
   plot( x, fx, 'LineWidth', 1.5), hold off;
   xlabel( "Puissance" );
   ylabel( "Densité de probabilité" );
   title( "PUISSANCE" ); 
   subtitle( "Densité de probabilité avec P_{bth}= " + string( round( 10*log10(Pbth_lin) ) ) + " dB" );  
   legend("Histogramme", "Densité théorique")
   grid on;                        


end
