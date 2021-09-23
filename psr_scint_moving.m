% The pulsar is moving at some velocity.
% 1. At each time, step, Gaussian generate a set of small angles that
% represent rays. These rays each pick up a phase and as it reaches Earth:
% 2. Add all the phases to give a single point (t,phase) at a constant k
% Repeat above until as set of points [t1,t2...,t] & [phase1,phase2...,]
% 3. Repeated the two steps for a range of K ----> DYNAMIC SPECTRUM
clear all;
clc;

psr_pos = [0,0]; screen_dist=2*3.086e16; observer=[4,0]*3.086e16; % [x,y]
velocity = [2e5,2e5] ; total_time = 1801; num_rays = 400;
%kmin=2*pi/(0.37474057*3e8); kmax=2*pi/(0.35269701*3e8);
kmin=2*pi/3e8; kmax=2.5*pi/3e8;

% psr_pos = [0,0]; screen_dist=50; observer=[100,0]; % [x,y]
% velocity = [0.001,0.001] ; total_time = 500; num_rays = 400;
% kmin=1e6; kmax=1.5e6;


% This is at each time
t_array = 1:total_time;
K_range = linspace(kmin,kmax,length(t_array));
K_array = 1:length(K_range);
phi_const = normrnd(0,pi/8,[1,length(num_rays)]);
%phi_const = linspace(-pi/8,pi/8,length(num_rays));
%phi_screen = normrnd(0,pi/1e8,[1,length(num_rays)]);
phase_array_at_k = zeros(length(K_array),length(t_array));
for K = K_array
    psr_pos=[0,0];
    for time_t = t_array
        psr_pos = psr_pos + velocity*(time_t-1);
        %phi = normrnd(0,pi/15e-7,[1,length(num_rays)]); % Here, ini ray ang. change every timestep as pulsar moves 
        phi = phi_const; % Here, the ini ray ang. is constant as pulsar moves
        r01 = (screen_dist - psr_pos(1)) ./ cos(phi);
        r12 = sqrt((r01.*sin(phi) + psr_pos(2))^2 + (observer(1)-screen_dist)^2);
        total_r = r01+r12;
        tot_phase_at_t = exp(1i*K.*total_r);
        phase_array_at_k(K,time_t) = sum(tot_phase_at_t);
        %tot_phase_at_t = exp(1i*K*sum(total_r));
        %phase_array_at_k(K,time_t) = tot_phase_at_t;
        
%         for rays = 1:num_rays
%             phi = normrnd(0,pi/8); % the angle of rays change every timestep as pulsar moves 
%             r01 = (screen_dist - psr_pos(1)) / cos(phi);
%             r12 = sqrt((r01*sin(phi) + psr_pos(2))^2 + (observer(1)-screen_dist)^2);
%             total_r = r01+r12;
%             %r01_phase = exp(1i*K*r01); r12_phase = exp(1i*K*r12);
%             tot_phase_at_t = exp(1i*K*total_r);
%             phase_array_temp(1,rays) = tot_phase_at_t;
%         end
        %phase_array_at_k(K,t) = sum(phase_array_temp); % interference of %all paths at time t, For loop
        
    end
    
end


%[X,Y] = meshgrid(t_array,K_range);
%surf(X,Y,real(phase_array_at_k))
figure(1)
surf(real(phase_array_at_k))
%shading interp
view(2)
% 
% figure(2)
% imshow(phase_array_at_k)

% sz=1;mkr='.';
% for ks = K_array
%     figure(2)
%     scatter(t_array,repelem((2*pi/(K_range(ks))),length(t_array)),sz,real(phase_array_at_k(ks,:)),mkr)
%     
%     shading interp
%     hold on
% end