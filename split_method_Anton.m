%Run cmap2.m first to set colormap

clear

%cmap

% load('Vdistor3.mat');
% load('z2.mat', 'z2')

z2 = double(rgb2gray(imread('abs_pot.jpg')));
% I2_ = rgb2gray(I_)
% z2 = double(I2_)


% Converting image to 2D-Matrix

V2 = double(imread('ISM_deformation_blur200.jpg'));
% V2 = double(rgb2gray(imread('ISM__.jpg')))

% I2 = rgb2gray(I)
% V2 = double(I2)

% VARIABLES

theta = 90

pulsing_period = 100

A=4e9;
V=A*V2/1000; % V is the potential
d=0.01;0.001; % damping parameter

NN=600; % number of steps

Nx=1024; % dimensions of the simulation
Ny=1024;

tdir=-1; % arrow of time (why is it negative?)

x0=0.0;
y0=-8.0; % Signal starts in middle top of simulation

x02=-8; % FUNCTION?
y02=0;

sprx=50;
spry=1;
xmin=-10;
xmax=10;10; % WHAT IS THE PURPOSE OF SECOND 10?
ymin=-10;-10;
ymax=10;10;
kk=10; % momentum

velocity = 1; % relative velocity of ISM cloud (in pixels per iteration)

grid=zeros(1, NN); % grid for primary spectrum

x_earth = 513; % position of earth (x-value)
y_earth = 800; % position of earth (y-value)

delx=(xmax-xmin)/(Nx-1); %  x-space grid spacing
dely=(ymax-ymin)/(Ny-1); %  y-space grid spacing

xi=xmin:delx:xmax; 
yi=ymin:dely:ymax;
[xt,yt] = meshgrid(xi,yi); % array giving the dimensions of the simulation 


% x1 is the initial wavefunction (Gaussian)

% x1= exp(-(xt-x0).^2/sprx-(yt-y0).^2/spry).*...
%    exp(.2*1i*kk*xt) ; % non-normalised 

% Initial wavefunction
x1=exp(-(xt-x0).^2/(sprx)-(yt-y0).^2/(spry)).*exp(1i*kk*(xt*cosd(theta)+yt*sind(theta))); % Gaussian
% x1=exp(-(xt-x0).^2/sprx-(yt-y0).^2/spry).*besselj(0,kk*sqrt((xt-x0).^2+(yt-y0).^2)); % Bessel beam
% Normalization of wavefunction
% norm1=sum(sum(conj(x1).*x1));
% x1=x1/sqrt(norm1);

norm1=sum(sum(conj(x1).*x1));
x1=x1/sqrt(norm1); % normalised



% momentum grid:
[pxt,pyt] = meshgrid( -pi/(xmax-xmin)*(Nx-1)*(Nx-1):2*pi/delx:pi/(xmax-xmin)*(Nx-1)*(Nx-1),  -pi/(ymax-ymin)*(Ny-1)*(Ny-1):2*pi/dely:pi/(ymax-ymin)*(Ny-1)*(Ny-1) ); % meshgrid is 'overlay' of two 2D-matrices

% the momentum propagation  matrix exp[-i p^2/2m dt]
tau1=-2e-9; % time increment for the evolution of the wave
avec=exp(-tau1*tdir*1i*( (pxt.^2+ pyt.^2))/2); 
% this takes care of some oddities of the way engineers do FFT:
avec(:,:) = circshift(avec(:,:),[round(Nx/2) round(Ny/2)]);
% the position propagation matrix exp[-i v(x) dt]
dvec=exp(-double(z2')*d)'; % absorption potential

% DEFINITION OF XVEC USED TO BE HERE!

vidObj = VideoWriter('PulsarVideo.avi');
vidObj.FrameRate=10;
vidObj.Quality=100;
open(vidObj);

x=x1;%As*x1;

t=0;
hbar=1;
w=1/abs(tau1);

E=1e8:1e8:5e8;%logspace(6,8,5);%0.5e6:1e8:5e8;%5e8;1e8;1e20;hbar*w; % energy, divided by hbar (correction: h not h bar) gives the frequency of the pulse
xeig=cell(1,length(E));

% creating a 1000x1000 array for each of the 5 entries of xeig
for i=1:length(E)
    xeig{i}=zeros(Ny,Nx);
end


for j=1:NN
    
if mod(j, pulsing_period) == 0   
    x = x + x1;
end    
    
    remaining=NN-j
    
    norm=sum(sum(conj(x).*x));
    x(:,:)  = x(:,:)/sqrt(norm);
    
    
    xvec=exp(-1i.*tdir*V'*tau1 )'.*dvec; % multiplication of propagators of absorption pot. and deformation pot.

    
    
    % THIS IS THE SPLIT-OPERATOR METHOD
    % propagate x part
    x(:,:)= xvec.*x(:,:);
    % FFT into momentum space
    x(:,:)= fft2(x(:,:));
    %propagate p part
    x(:,:) = avec.*x(:,:);
    % back to coordinate space
    x(:,:)= ifft2(x(:,:));
    
    % propagating the eigenvalues at all five energies
    for i=1:length(E)
        xeig{i}=xeig{i}+x.*exp(-1i.*(j-1)*tdir*E(i)'*tau1 )';
        
    end
    
    figure(2)
    imagesc(real(x))
    axis equal off
    colormap(uu)
    top=max(max(abs(real(x))));
    caxis([-top top])
    sq(j)=getframe;
    writeVideo(vidObj,sq(j));
    
    
    
    % movement of ISM cloud (move the potential towards the right)
    
    % xvec = circshift(xvec, 1, 2)
    V = circshift(V, velocity, 2);
    
    % PRIMARY SPECTRUM (measurement of intensity)
    
    
    grid(1, j)=x(y_earth, x_earth);
    
    
end

close(vidObj);




figure(),pcolor(xi,yi,V'),shading flat, colormap(flipud(gray)),axis equal tight

for i=1:length(E)
    figure();
    imagesc(xi,yi,real(xeig{i}))
    top=max(max(abs(real(xeig{i}))));
    caxis([-top top])
    shading flat, axis equal tight,colormap(flipud(uu))
end