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


t_array = linspace(1,total_time,total_time); % Time in seconds, create array of time
K_range = linspace(kmin,kmax,length(t_array)); % the Range of wavewector k that pulsar emits
K_array = 1:length(K_range); % 'Index' for K_range to store k values into arrays
phi_const = normrnd(0,pi/8,[1,length(num_rays)]); % initialize a set of angles from which the ray is sent

phase_array_at_k = zeros(length(K_array),length(t_array)); % This matrix is used to store phase
for K = K_array
    psr_pos=[0,0]; %The position of pulsar is 'reset' for each k.
    
    for time_t = t_array
    
        psr_pos = psr_pos + velocity*(time_t-1); %update position of pulsar at every timestep
        
        phi = phi_const; % Retrieve angle from intially generated angles
        r01 = (screen_dist - psr_pos(1)) ./ cos(phi); % Distance from pulsar to screen
        r12 = sqrt((r01.*sin(phi) + psr_pos(2))^2 + (observer(1)-screen_dist)^2); % Distance screen to observer
        total_r = r01+r12; % Total distance
        tot_phase_at_t = exp(1i*K_range(K).*total_r); % Phase of the distance at this specific timestep
        phase_array_at_k(K,time_t) = sum(tot_phase_at_t); % Add up all the phases at this time specific timestep: Interference

        
    end
    
end



figure(1)
surf(real(phase_array_at_k))
shading interp
view(2)
