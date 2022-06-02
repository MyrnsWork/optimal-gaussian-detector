% Maël Laviec
% mlaviec@enseirb-matmeca.fr
% ENSEIRB-MATMECA
% Juin, 2022

function [ logLRT_lin,...
           detectionMap,...
           gammaLogLRT_dB  ] = optimalDetector( imagetteChannelIQ_lin,...
                                                R,...
                                                Pfa,...
                                                targetIQ,...
                                                steringVector,...
                                                typeTarget,...
                                                dimensionFreq           )
                               
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
