%% Initialisation
clc; 
close all;
clear;
dbstop if error;


%% Hyperparamètres
kB = 1.38e-23;
c  = 3e8;

%% Paramètres
%géométrie des cartes
Nrec  = 128;                 %nombre de récurrences, scalaire
Ncd   = 1000;                %nombre de cases distances, scalaire
Ncell = Nrec * Ncd;          %nombre de cellule, scalaires

%radar
fe            = 10e9;        %fréquence de la porteuse, scalaire [Hz]
lambda_m      = c / fe;      %longueur d'onde de la porteuse, scalaire [m]
Bp            = 200e6;       %bande passante, scalaire [m]
Li            = 1/Bp;        %longueur d'impulsion, scalaire [s]  
theta_3dB_deg = 15;
theta_3dB_rad = theta_3dB_deg * pi/180;
wrot_rads     = 20 * pi/30;  %vitesse de rotation du radar, scalaire [rad/s] 
PRF           = 1500;        %fréquence de récurrence, scalaire [Hz]

%caractérisation du bruit thermique
T_degCel = 20;               %température en degrés celsus
T_K      = 273.15 + T_degCel;%température en degrés Kelvin
F_dB     = 6;
F_lin    = 10^( F_dB/10);    %facteur de bruit de l'électronique analogique
Pbth_lin = kB * T_K * Bp * F_lin;
                             %Puissance du bruit thermique
Pbth_dB  = 10*log10( Pbth_lin );                              

%caractéristiques de la cible
SNR_dB          = 20;
speedTarget     = 5;
targetFrequency = 2 * speedTarget / lambda_m;
swerlingType    = 1;
azimut0_deg     = 90;
azimut0_rad     = azimut0_deg * pi/180;


%% Imagette de bruit blanc complexe (channels I et Q)
[ imagetteAmplitude_lin,...
  ~                        ] = createImagette( Pbth_lin,...
                                               Ncd,...
                                               Nrec        );
                                           
puissanceImagette_dB = 10*log10( 1/Ncell * sum( abs(imagetteAmplitude_lin(:) ).^2 ) );

%% Ajout de la cible 
[ imagetteAmplitude_lin,...
  imagettePuissance_lin,...
  rangeGate                ] = addTarget( imagetteAmplitude_lin,...
                                          SNR_dB,... 
                                          swerlingType,...
                                          targetFrequency,...
                                          azimut0_rad,...
                                          wrot_rads,...
                                          theta_3dB_rad,...
                                          PRF,...
                                          Bp,...
                                          Li,...
                                          Pbth_dB,...
                                          Ncd,...
                                          Nrec                     );

%% Carte Doppler de l'imagette
timeWindow = repmat( hamming(Nrec, 'periodic')', Ncd, 1);
carteDopplerAmplitude_lin = 1/Nrec * fftshift( fft(imagetteAmplitude_lin .* timeWindow, [], 2), 2 );
carteDopplerPuissance_lin = abs(carteDopplerAmplitude_lin).^2;

puissanceCarteDoppler_dB = 10*log10( 1/Ncd * sum( abs(carteDopplerAmplitude_lin(:)).^2 ) );

%% Figures
fig1 = displayImagette( imagettePuissance_lin,...
                        1                        );
fig2 = displayCarteDoppler( carteDopplerPuissance_lin,...
                            2                            );
                        

%% Fonctions 
function [ imagetteAmplitude_lin,...
           imagettePuissance_lin    ] = createImagette( Pbth_lin,...
                                                        Ncd,...
                                                        Nrec        )
                                                    
    imagetteAmplitude_lin = sqrt( Pbth_lin/2 ) * ( randn(Ncd, Nrec) + 1j*randn(Ncd, Nrec) );
    imagettePuissance_lin = abs(imagetteAmplitude_lin).^2;
                                                    
end

function fig = displayImagette( imagettePuissance_lin,...
                                iFigure                  )

    fig = figure(iFigure);
    imagesc( 10*log10( abs(imagettePuissance_lin) ) );
    colorbar, colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Réccurence')
    ylabel('Distance relative')
    title('Imagette de puissance [dB]')


end

function fig = displayCarteDoppler( carteDopplerPuissance_lin,...
                                    iFigure                      )

    fig = figure(iFigure);
    imagesc( 10*log10( abs(carteDopplerPuissance_lin) ));
    colorbar, colormap('parula')
    set(gca, 'YDir', 'normal')
    xlabel('Fréquence')
    ylabel('Distance relative')
    title('Carte Doppler de puissance [dB]')


end



function [ imagetteAmplitude_lin,...
           imagettePuissance_lin,...
           rangeGate                ] = addTarget( imagetteAmplitude_lin,...
                                                   SNR_dB,... 
                                                   swerlingType,...
                                                   targetFrequency,...
                                                   azimut0_rad,...
                                                   wrot_rads,...
                                                   theta_3dB_rad,...
                                                   Fr,...
                                                   Bp,...
                                                   Li,...
                                                   Pbth_dB,...
                                                   Ncd,...
                                                   Nrec                     )
                                            

    c = 3e8;
    rtot = 10;
    
    Tr = 1/Fr;  
    t  = ((1 : Nrec) - Nrec/2) * Tr;
    
    rangeGateTot = floor((rtot * (2*Bp) / c) / 2);
    rangeGateTot = (-rangeGateTot: rangeGateTot)'; 
    rangeGate = randi( [ ceil(length(rangeGateTot)/2),...
                         Ncd - ceil(length(rangeGateTot)/2) ] );
    rangeGateTot = rangeGateTot + rangeGate;                
    azimut_rad = t.* wrot_rads + azimut0_rad;
    r = (rangeGateTot - rangeGate) * c/(2*Bp);
    
    Baz = exp( - (azimut_rad - azimut0_rad).^2 ./ ( 2 * 2*log(2)*theta_3dB_rad.^2 ) );
    %% Parameters 
    mu = Bp / (2 * Li);
    
    %% Vectors (fastTime : range dimension,  slowTime : azimut dimension
    fastTime = 2 * r / c/4;
    slowTime = targetFrequency;
    
    %% Ambiguity function
    val1   = 1. - abs(fastTime) / Li;
    val2   = Li .* val1;
    val3   = (slowTime + mu .* fastTime);
    val4   = val2 .* val3;
    Brange = abs( val1 .* sinc( val4 ) ).^2;

    switch swerlingType
        case 0
            swerling = 10.^( (Pbth_dB + SNR_dB)/10 );
        case 1
            swerling = exprnd( 10.^( (Pbth_dB + SNR_dB)/10 ), length(rangeGateTot), Nrec );
            
    end
    target   = sqrt( swerling .* Baz .* Brange) .* exp(1j*2*pi*targetFrequency * t);
    imagetteAmplitude_lin(rangeGateTot, :) = imagetteAmplitude_lin(rangeGateTot, :) + target;
    
    imagettePuissance_lin = abs(imagetteAmplitude_lin).^2;
    
end



