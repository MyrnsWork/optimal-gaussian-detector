% Maël Laviec
% mlaviec@enseirb-matmeca.fr
% ENSEIRB-MATMECA
% Juin, 2022


%% Initialisation
clc; 
close all;
clear;
dbstop if error;


%% Options du script
displayFigures   = true;                                                   % affichage des figures, 1x1
isCovarianceKown = false;                                                  % connaissance de la matrice de covariance du bruit thermique, 1x1


%% Hyperparamètres
kB  = 1.38e-23;                                                            % constante de Boltzmann, 1x1 [m2 kg s-2 K-1]
c   = 3e8;                                                                 % célérité de la lumière, 1x1 [m/s]
Pfa = 1e-4;                                                                % probabilité de fausse alarme, 1x1


%% Paramètres
% géométrie des cartes
Nrec  = 128;                                                               % nombre de récurrences, 1x1
Ncd   = 150;                                                               % nombre de cases distance, 1x1
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
SNR_dB          = 9;                                                       % rapport signal sur bruit, 1x1                                                                          
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



%% Calcul de la matrice de covariance
scmR = calcCovarianceMatrix( imagetteChannelIQ_lin,...
                             R,...
                             isCovarianceKown,...
                             2                        );


%% Détecteur optimal
[ logLRT_lin,...
  detectionMap,...
  gammaLogLRT_dB  ] = optimalDetector( imagetteChannelIQWithTarget_lin,...
                                       scmR,...
                                       Pfa,...
                                       targetIQ,...
                                       steringVector,...
                                       typeTarget,...
                                       1                                  ) ;    



%% Figures
if displayFigures
    fig1 = displayImagettes( imagettePuissance_lin,...
                             imagettePuissanceWithTarget_lin,...
                             detectionMap,...
                             rangeIndex,...
                             Nrec,...
                             1                                  );     
    
    
    fig2 = displayHistogram( imagetteAmplitude_lin,...
                             imagettePuissance_lin,...
                             Pbth_lin,...
                             2                        );
                         
    
    
    fig3 = displayLogLRT( logLRT_lin,...
                          gammaLogLRT_dB,...
                          Ncd,...
                          3                 );                 

end
                      


%% Fonctions
% calculs
function [ imagetteChannelIQ_lin,...
           imagetteAmplitude_lin,...
           imagettePuissance_lin,...
           rangeIndex               ] = addTarget( imagetteChannelIQ_lin,...
                                                   targetIQ,...
                                                   dimensionFreq            )
                                                         
   switch dimensionFreq
       case 1
           rangeIndex = randi( size(imagetteChannelIQ_lin, 1) );
           imagetteChannelIQ_lin(:, rangeIndex) = imagetteChannelIQ_lin(:, rangeIndex) + targetIQ;  
          
       case 2
           rangeIndex = randi( size(imagetteChannelIQ_lin, 2) );
           imagetteChannelIQ_lin(rangeIndex, :) = imagetteChannelIQ_lin(rangeIndex, :) + targetIQ;  
           
       otherwise
           error( "Dimension inconnue ou non supportée" );
           
   end
   
   
   imagetteAmplitude_lin = abs( imagetteChannelIQ_lin );    %non optimisé, il suffit de faire sur seulement la ligne concernée
   imagettePuissance_lin = abs( imagetteChannelIQ_lin ).^2; %non optimisé, il suffit de faire sur seulement la ligne concernée
   
end

function R = calcCovarianceMatrix( imagetteChannelIQ_lin,...
                                   R,...
                                   isCovarianceKown,...
                                   dimensionRec             )

    if ~isCovarianceKown                                                     
       Ncd = size(imagetteChannelIQ_lin, dimensionRec);
       R   = 1/Ncd * pagemtimes(imagetteChannelIQ_lin, imagetteChannelIQ_lin');
       R   = diag( diag( R ) );

    end
        

end

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
            Tr              = 1/PRF;
            amplitude       = sqrt(SNR_lin * Pbth_lin);
            phase           = 0;
            targetAmplitude = amplitude .* exp(1j .* phase);
            t               = (0 : Tr : (Nrec-1) * Tr)';  
            steringVector   = exp( 1j*2*pi * targetFrequency * t);
            targetIQ        = targetAmplitude .* steringVector; 

        case "swerling1"                                                   % amplitude suivant une loi de Rayleigh, phase suivant une loi uniforme
            Tr              = 1/PRF;
            targetPower     = SNR_lin * Pbth_lin;
            sigmaTarget     = sqrt( targetPower/2 );
            targetAmplitude = sigmaTarget * ( randn + 1j*randn );          % hypothèse : le tirage d'ampltitude est 
                                                                           % constant sur le temps d'intégration (Tint = Nrec/PRF )
            t               = (0 : Tr : (Nrec-1) * Tr)';  
            steringVector   = exp( 1j*2*pi * targetFrequency * t);
            targetIQ        = targetAmplitude .* steringVector;   

        case {"swerling0", "swerling5"}                                    % amplitude constante, phase suivant une loi uniforme
            Tr              = 1/PRF;
            targetPower     = SNR_lin * Pbth_lin;
            Phi             = 2 * pi * rand;                               % phase aléatoire uniforme sur l'intervalle [0;2pi]
            sigmaTarget     = sqrt( targetPower/2 );
            targetAmplitude = sigmaTarget * exp(1j*Phi);                   % hypothèse : le tirage d'ampltitude est 
                                                                           % constant sur le temps d'intégration (Tint = Nrec/PRF )
            t               = (0 : Tr : (Nrec-1) * Tr)'; 
            steringVector   = exp( 1j*2*pi * targetFrequency * t);
            targetIQ        = targetAmplitude .* steringVector;   

       case "unknow"
            Tr              = 1/PRF;
            SNR_dB          = rand * 30 - 5;
            SNR_lin         = 10^( SNR_dB/10 );
            targetPower     = SNR_lin * Pbth_lin;                          % amplitude déterministe et inconnue
            sigmaTarget     = sqrt( targetPower/2 );
            targetAmplitude = sigmaTarget * ( randn + 1j*randn );          % hypothèse : le tirage d'ampltitude est 
                                                                           % constant sur le temps d'intégration (Tint = Nrec/PRF )
            t               = (0 : Tr : (Nrec-1) * Tr)';   
            steringVector   = exp( 1j*2*pi * targetFrequency * t);
            targetIQ        = targetAmplitude .* steringVector;               

        otherwise
            error( "type de cible inconnu ou non supporté" );
    end

                               
                               
                               
end

function [ logLRT_lin,...
           detectionMap,...
           gammaLogLRT_dB  ] = optimalDetector( imagetteChannelIQ_lin,...
                                                R,...
                                                Pfa,...
                                                targetIQ,...
                                                steringVector,...
                                                typeTarget,...
                                                dimensionFreq            )
                               
   switch typeTarget

     case "deterministic"
        logLRT_lin      = 2 * real( pagemtimes(targetIQ'/R, imagetteChannelIQ_lin) );
        sigma0          = sqrt( real( 2 * targetIQ' * (R \ targetIQ) ) );
        gammaLogLRT_lin = sigma0 * qfuncinv(Pfa);
        gammaLogLRT_dB  = 10 * log10(gammaLogLRT_lin);  
        
     case {"swerling1", "swerling0", "swerling5", "unknow"}  
        logLRT_lin      = abs( pagemtimes(steringVector'/R, imagetteChannelIQ_lin) ).^2;
        sigma0          = real( steringVector' * (R \ steringVector) );                           % partie réelle uniquement (à cause des erreurs numériques)
        gammaLogLRT_lin = -sigma0 * log( Pfa );
        gammaLogLRT_dB  = 10 * log10(gammaLogLRT_lin);  

    otherwise
        error( "type de cible inconnu ou non supporté" );

   end

   detectionMap = (logLRT_lin >= gammaLogLRT_lin);

   switch dimensionFreq
       case 1
           Nrec         = size(imagetteChannelIQ_lin, 1); 
           detectionMap = repmat(detectionMap, Nrec, 1);
          
       case 2
           Nrec         = size(imagetteChannelIQ_lin, 2); 
           detectionMap = repmat(detectionMap, 1, Nrec);
           
       otherwise
           error( "Dimension inconnue ou non supportée" );
           
   end

end


% affichage
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
    xlabel('Case distance')
    title('Imagette de puissance sans cible [dB]')
    
    %----------------------------------------------------------------------
    subplot 132
    imagesc( 10*log10( abs(imagettePuissanceWithTarget_lin') ), [minPower maxPower]);
    colorbar, colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Réccurence')
    xlabel('Case distance')
    title('Imagette de puissance avec cible [dB]')

    %----------------------------------------------------------------------
    subplot 133
    imagesc( detectionMap' ), hold on;
    L = line(ones(2), ones(2), 'LineWidth',2);      
    cmap = parula(2);  
    set(L,{'color'}, mat2cell(cmap,ones(1,2),3))    
    plot( (0 : Nrec+1), rangeTarget * ones(1, Nrec+2), 'r', 'LineWidth', 1.5 )
    hold off           
    legend( 'Non-détection', 'Détection', 'Index de la cible' )    
    colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Réccurence')
    xlabel('Case distance')
    title('Carte de détection')

end      

function fig = displayLogLRT( logLRT_lin,...
                              gammaLogLRT_dB,...
                              Ncd,...
                              iFigure           )

    fig = figure( iFigure );
    plot( 10*log10( abs( logLRT_lin(1,:) ) ), 'LineWidth', 1.5 ), hold on;
    plot( gammaLogLRT_dB * ones(1, Ncd), 'LineWidth', 1.5 ), hold off;
    legend( "LogLRT", "Seuil Optimal" )
    ylabel('Rapport du log-vraisemblance [dB]')
    xlabel('Case distance')
    xlim([1 Inf])
    title('Evolution du rapport de log-vraisemblance en fonction de la distance [dB]')
    set(fig, 'Units', 'Normalized', 'Position', [0 0 1 1]);


end      







