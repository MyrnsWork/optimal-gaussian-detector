function fig = displayLogLRT( logLRT_lin,...
                              gammaLogLRT_dB,...
                              Ncd,...
                              iFigure           )

    fig = figure( iFigure );
    plot( 10*log10( abs( logLRT_lin(1,:) ) ), 'LineWidth', 1.5 ), hold on;
    plot( gammaLogLRT_dB * ones(1, Ncd), 'LineWidth', 1.5 ), hold off;
    legend( "LogLRT", "Seuil Optimal" )
    ylabel('Rapport du log-vraisemblance [dB]')
    xlabel('Distance relative')
    xlim([1 Inf])
    title('Evolution du rapport de log-vraisemblance [dB]')
    set(fig, 'Units', 'Normalized', 'Position', [0 0 1 1]);


end      
