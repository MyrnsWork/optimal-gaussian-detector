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