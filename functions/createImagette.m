function [ imagetteChannelIQ_lin,...
           imagetteAmplitude_lin,...
           imagettePuissance_lin    ] = createImagette( Pbth_lin,...
                                                        Ncd,...
                                                        Nrec        )
                                                    
    imagetteChannelIQ_lin = sqrt( Pbth_lin/2 ) * ( randn(Nrec, Ncd) + 1j*randn(Nrec, Ncd) );
    imagetteAmplitude_lin = abs(imagetteChannelIQ_lin);
    imagettePuissance_lin = abs(imagetteChannelIQ_lin).^2;
                                                    
end