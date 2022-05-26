function [ ambiguityFunction,...
           fastTime,...
           slowTime             ] = lfm_ambg( Li,... 
                                              Bp,...
                                              N,...
                                              M     )
    %% Parameters 
    mu = Bp / (2 * Li);
    
    %% Vectors (fastTime : range dimension,  slowTime : azimut dimension
    fastTime = (-Li : 2*Li/N : Li);
    slowTime = (-Bp : 2*Bp/M : Bp)';
    
    %% Ambiguity function
    val1 = 1. - abs(fastTime) / Li;
    val2 = Li .* val1;
    val3 = (slowTime + mu .* fastTime);
    val4 = val2 .* val3;
    ambiguityFunction = abs( val1 .* sinc( val4 ) ).^2;

end
