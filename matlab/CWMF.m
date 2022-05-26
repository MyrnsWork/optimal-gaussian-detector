%% Initialisation
clc; 
close all;
clear;
dbstop if error;

%% Hyperparamètres
kB  = 1.38e-23;                                                            % constante de Boltzmann, 1x1 [m2 kg s-2 K-1]
c   = 3e8;                                                                 % célérité de la lumière, 1x1 [m/s]
Pfa = 1e-4;                                                                % probabilité de fausse alarme, 1x1

%% Paramètres
% géométrie des cartes
Nrec  = 64;                                                                % nombre de récurrences, 1x1
Ncd   = 100;                                                               % nombre de cases distance, 1x1
Ncell = Nrec * Ncd;                                                        % nombre de cellules, 1x1

% radar
fe     = 10e9;                                                             % fréquence de la porteuse, 1x1 [Hz]
lambda = c/fe;                                                             % longueur d'onde de la porteuse, 1x1 [m]
Bp     = 200e6;                                                            % bande passante, 1x1 [m]
Li     = 1/Bp;                                                             % longueur d'impulsion, 1x1 [s]  
PRF    = 2000;                                                             % fréquence de récurrence, 1x1 [Hz]

% caractérisation du bruit thermique
T_degCel = 20;                                                             % température en degrés Celsus, 1x1 [°C]
T_K      = 273.15 + T_degCel;                                              % température en degrés Kelvin, 1x1 [°K]
F_dB     = 6;                                                              % facteur de bruit de l'électronique analogique, 1x1
F_lin    = 10^( F_dB/10 );                                                 
Pbth_lin = kB * T_K * Bp * F_lin;                                          % puissance du bruit thermique, 1x1 [W]                                                                                                   
Pbth_dB  = 10*log10( Pbth_lin );                                           
R        = Pbth_lin * eye(Nrec);                                           % matrice de covariance du bruit thermique, Nrec x Nrec

% caractéristiques de la cible
SNR_dB          = 10;                                                      % rapport signal sur bruit, 1x1                                                                          
SNR_lin         = 10^( SNR_dB/10 );
typeTarget      = "swerling1";                                             % type de fluctuations de la cible, 1x1   
speedTarget     = 5;                                                       % vitesse radiale de la cible, 1x1 [m/s]   
targetFrequency = 2 * speedTarget / lambda;                                % fréquence Doppler de la cible, 1x1 [Hz]   


%% Imagette de bruit blanc complexe (channels I et Q)
[ imagetteChannelIQ_lin,...
  imagetteAmplitude_lin,...
  imagettePuissance_lin    ] = createImagette( Pbth_lin,...
                                               Ncd,...
                                               Nrec        );
                                           
%% Création de la cible 
[ targetIQ,...
  targetAmplitude,...
  steringVector      ] = createTarget( SNR_lin,...
                                       Pbth_lin,...
                                       targetFrequency,...
                                       PRF,...
                                       typeTarget,...
                                       Nrec               );

%% Ajout de la cible
[ imagetteChannelIQWithTarget_lin,...
  imagetteAmplitudeWithTarget_lin,...
  imagettePuissanceWithTarget_lin,...
  rangeIndex                         ] = addTarget( imagetteChannelIQ_lin,...
                                                    targetIQ,...
                                                    1                        );                            

%% Détecteur optimal
[ logLRT_lin,...
  detectionMap,...
  gammaLogLRT_dB  ] = optimalDetector( imagetteChannelIQWithTarget_lin,...
                                       R,...
                                       Pfa,...
                                       targetIQ,...
                                       steringVector,...
                                       typeTarget,...
                                       1                                  );    


%% Figures
fig1 = displayImagettes( imagettePuissance_lin,...
                         imagettePuissanceWithTarget_lin,...
                         detectionMap,...
                         rangeIndex,...
                         Nrec,...
                         1                        );     


fig2 = displayHistogram( imagetteAmplitude_lin,...
                         imagettePuissance_lin,...
                         Pbth_lin,...
                         2                        );
                     


fig3 = displayLogLRT( logLRT_lin,...
                      gammaLogLRT_dB,...
                      Ncd,...
                      3                 );                 


                      


%% Fonctions 
function [ imagetteChannelIQ_lin,...
           imagetteAmplitude_lin,...
           imagettePuissance_lin    ] = createImagette( Pbth_lin,...
                                                        Ncd,...
                                                        Nrec        )
                                                    
    imagetteChannelIQ_lin = sqrt( Pbth_lin/2 ) * ( randn(Nrec, Ncd) + 1j*randn(Nrec, Ncd) );
    imagetteAmplitude_lin = abs(imagetteChannelIQ_lin);
    imagettePuissance_lin = abs(imagetteChannelIQ_lin).^2;
                                                    
end

function [ targetIQ,...
           targetAmplitude,...
           steringVector      ] = createTarget( SNR_lin,...
                                                Pbth_lin,...
                                                targetFrequency,...
                                                PRF,...
                                                type,...
                                                Nrec               )
                               
                               
    switch type
        case "deterministic"
            amplitude       = sqrt(SNR_lin * Pbth_lin);
            phase           = 0;
            targetAmplitude = amplitude .* exp(1j .* phase);
            t               = (0 : Nrec-1)' * 1/PRF;  
            steringVector   = exp( 1j*2*pi * targetFrequency * t);
            targetIQ        = targetAmplitude .* steringVector; 

        case "swerling1"
            targetPower     = SNR_lin * Pbth_lin;
            sigmaTarget     = sqrt( targetPower/2 );
            targetAmplitude = sigmaTarget * ( randn + 1j*randn );          % hypothèse : le tirage d'ampltitude est 
                                                                           % constant sur le temps d'intégration (Tint =9 Nrec/PRF )
            t               = (0 : Nrec-1)' * 1/PRF;  
            steringVector   = exp( 1j*2*pi * targetFrequency * t);
            targetIQ        = targetAmplitude .* steringVector;   

        otherwise
            errror( "type de cible inconnu ou non supporté" );
    end

                               
                               
                               
end

function [ imagetteChannelIQ_lin,...
           imagetteAmplitude_lin,...
           imagettePuissance_lin,...
           rangeIndex               ] = addTarget( imagetteChannelIQ_lin,...
                                                   targetIQ,...
                                                   dimensionRec             )
                                                         
   if dimensionRec == 1                                                      
       rangeIndex = randi( size(imagetteChannelIQ_lin, 2) );
       imagetteChannelIQ_lin(:, rangeIndex) = imagetteChannelIQ_lin(:, rangeIndex) + targetIQ;  
      
   elseif dimensionRec == 2
       rangeIndex = randi( size(imagetteChannelIQ_lin, 1) );
       imagetteChannelIQ_lin(rangeIndex, :) = imagetteChannelIQ_lin(rangeIndex, :) + targetIQ;  
       
   else 
       error( "Dimension inconnue ou non supportée" );
       
   end
   
   
   imagetteAmplitude_lin = abs( imagetteChannelIQ_lin );    %non optimisé, il suffit de faire sur seulement la ligne concernée
   imagettePuissance_lin = abs( imagetteChannelIQ_lin ).^2; %non optimisé, il suffit de faire sur seulement la ligne concernée
   
end

function [ logLRT_lin,...
           detectionMap,...
           gammaLogLRT_dB  ] = optimalDetector( imagetteChannelIQ_lin,...
                                                R,...
                                                Pfa,...
                                                targetIQ,...
                                                steringVector,...
                                                typeTarget,...
                                                dimensionRec             )
                               
   if dimensionRec == 1                                                      
       Ncd = size(imagetteChannelIQ_lin, 2);
%        Nrec = size(imagetteChannelIQ_lin, 1); 
      
   elseif dimensionRec == 2
       Ncd = size(imagetteChannelIQ_lin, 1);
%        Nrec = size(imagetteChannelIQ_lin, 2); 
       
   else 
       error( "Dimension inconnue ou non supportée" );
       
   end
   
   logLRT_lin = zeros( size(imagetteChannelIQ_lin) );
   for iRangeGate = (1  : Ncd)
      switch typeTarget

        case "deterministic"
            logLRT_lin(:, iRangeGate) = 2 * real(targetIQ' * (R \ imagetteChannelIQ_lin(:, iRangeGate)) );
            sigma0                    = sqrt( real( 2 * targetIQ' * (R \ targetIQ) ) );
            gammaLogLRT_lin           = sigma0 * qfuncinv(Pfa);
            gammaLogLRT_dB            = 10 * log10(gammaLogLRT_lin);  
            
        case "swerling1"
            logLRT_lin(:, iRangeGate) = abs(steringVector' * (R \ imagetteChannelIQ_lin(:, iRangeGate)) ).^2;
            sigma0                    = real( steringVector' * (R \ steringVector) );                           % partie réelle uniquement (à cause des erreurs numériques)
            gammaLogLRT_lin           = -sigma0 * log( Pfa );
            gammaLogLRT_dB            = 10 * log10(gammaLogLRT_lin);  

        otherwise
            errror( "type de cible inconnu ou non supporté" );

      end

   end

    detectionMap = logLRT_lin >= gammaLogLRT_lin;

end

function fig = displayImagettes( imagettePuissance_lin,...
                                 imagettePuissanceWithTarget_lin,...
                                 detectionMap,...
                                 rangeTarget,...
                                 Nrec,...
                                 iFigure                  )

    fig = figure(iFigure);
    set(fig, 'Units', 'Normalized', 'Position', [0 0 1 1]);
    minPower = 10 * log10( min( min( imagettePuissance_lin(:) ), min( imagettePuissanceWithTarget_lin(:) ) ) );
    maxPower = 10 * log10( max( max( imagettePuissance_lin(:) ), max( imagettePuissanceWithTarget_lin(:) ) ) );

    %----------------------------------------------------------------------
    subplot 131
    imagesc( 10*log10( abs(imagettePuissance_lin') ), [minPower maxPower]);
    colorbar, colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Réccurence')
    ylabel('Distance relative')
    title('Imagette de puissance sans cible [dB]')
    
    %----------------------------------------------------------------------
    subplot 132
    imagesc( 10*log10( abs(imagettePuissanceWithTarget_lin') ), [minPower maxPower]);
    colorbar, colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Réccurence')
    ylabel('Distance relative')
    title('Imagette de puissance avec cible [dB]')

    %----------------------------------------------------------------------
    subplot 133
    imagesc( detectionMap' ), hold on;
    plot( (1 : Nrec), rangeTarget * ones(1, Nrec), 'ro' ), hold off; 
    colorbar, colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Réccurence')
    ylabel('Distance relative')
    title('Carte de détection')

end      

function fig = displayLogLRT( logLRT_lin,...
                              gammaLogLRT_dB,...
                              Ncd,...
                              iFigure           )

    fig = figure(iFigure);
    plot( 10*log10( abs( logLRT_lin(1,:) ) ), 'LineWidth', 1.5 ), hold on;
    plot( gammaLogLRT_dB * ones(1, Ncd), 'LineWidth', 1.5 ), hold off;
    legend( "LogLRT", "Seuil Optimal" )
    ylabel('Rapport du log-vraisemblance [dB]')
    xlabel('Distance relative')
    title('Evolution du rapport de log-vraisemblance [dB]')
    set(fig, 'Units', 'Normalized', 'Position', [0 0 1 1]);


end      

function fig = displayHistogram( imagetteAmplitude_lin,...
                                 imagettePuissance_lin,...
                                 Pbth_lin,...
                                 iFigure                  )
                             
                             
                             
   fig = figure( iFigure );
   set(fig, 'Units', 'Normalized', 'Position', [0 0 1 1]);
   subplot 121
   N    = 1000;
   Pmax = 1 - 1/length( imagetteAmplitude_lin(:) );
   xMax = raylinv( Pmax, sqrt(Pbth_lin/2) );
   x    = (0 : N-1) * xMax/N;                                  
   fx   = raylpdf( x, sqrt(Pbth_lin/2) );                                   
                                     
   histogram( imagetteAmplitude_lin(:), 'Normalization', 'pdf'), hold on;
   plot( x, fx, 'LineWidth', 1.5), hold off;
   xlabel( "Amplitude [V]" );
   ylabel( "Densité de probabilité" );
   title( "Densité de probabilité de l'imagette pour P_{bth}=" + string( round( 10*log10(Pbth_lin) ) ) + "dB" ); 
   grid on;
   
   subplot 122
   N    = 1000; 
   Pmax = 1 - 1/length( imagettePuissance_lin(:) );
   xMax =  -Pbth_lin * log( 1 - Pmax ) ;
   x    = (0 : N-1) * xMax/N;                                  
   fx   = 1/Pbth_lin * exp( -x/Pbth_lin );                                   
                                     
   histogram( imagettePuissance_lin(:), 'Normalization', 'pdf'), hold on;
   plot( x, fx, 'LineWidth', 1.5), hold off;
   xlabel( "Puissance [W]" );
   ylabel( "Densité de probabilité" );
   title( "Densité de probabilité de l'imagette pour P_{bth}=" + string( round( 10*log10(Pbth_lin) ) ) + "dB" ); 
   grid on;                             
                                
                                
end


