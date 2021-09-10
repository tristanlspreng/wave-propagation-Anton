% learning

%% matrix = 1:1:5;



% xvec_abs = real(xvec*100+100)

% image(xvec_abs)

% avec_abs = imag(avec*100+100+100i)
% image(avec_abs)


%% AA = [1 4 7 10
%%      2 5 8 11
%%      3 6 9 12]
  
%% index = [4 1 2 3]

%% BB = circshift(AA, 1, 2)


% propagation of the 1000x1000 matrix (i.e. the potential)

% prop = zeros(1000, 1000)


  
 % Attempt to add circles around pulsar and earth
 
 % hold on

 % plot(1.0, 0.5, 'bo', 'MarkerSize', 50)
 
 % legend('sine', 'circle')
  
 % hold off
 
 
 % image(real(x*100000)) % xvec is not the wavefunction, but x is!
 
 
 plot(abs(grid));
 % title('Real part of wavefunction receaved at earth vs. time')
 % xlabel('No. of iterations')
 % ylabel('Im({\Psi})')
 
 % GET SECONDARY SPECTRUM
 
 figure;spectrogram(grid,128,120,128,1e3, 'yaxis');
 
 % figure;spectrogram(grid, 'yaxis');
 
 
 % with Blackman windowing function:
 % spectrogram(grid,blackman(128),120,128,1e3, 'yaxis') 

 
 % view(-45,65)
 colormap bone;
 
 % image(z2)
 

 

 
 
