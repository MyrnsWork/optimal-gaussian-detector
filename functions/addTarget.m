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
