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
