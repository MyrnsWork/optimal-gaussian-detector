clear
Bp = 200e6;
Li = 1/Bp;
N = 400;
M = N * 10;
[ ambiguityFunction,...
  fastTime,...
  slowTime             ] = lfm_ambg( Li,... 
                                     Bp,...
                                     N,...
                                     M     );
uncertaintyFunction = sqrt(ambiguityFunction);

figure(1)
mesh(fastTime, slowTime, ambiguityFunction)
xlabel ('Delay - seconds')
ylabel ('Doppler - Hz')
zlabel ('Ambiguity function')

figure(2)
mesh(fastTime, slowTime, ambiguityFunction)
xlabel ('Delay - seconds')
ylabel ('Doppler - Hz')

figure(3)
mesh(fastTime, slowTime, uncertaintyFunction)
xlabel ('Delay - seconds')
ylabel ('Doppler - Hz')
zlabel ('Uncertainty function')

figure(4)
contour(fastTime, slowTime,uncertaintyFunction)
xlabel ('Delay - seconds')
ylabel ('Doppler - Hz')

