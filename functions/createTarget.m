% Maël Laviec
% mlaviec@enseirb-matmeca.fr
% ENSEIRB-MATMECA
% Juin, 2022

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
