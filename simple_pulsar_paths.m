clear all;
clc;
% 

source = 0; screen_dist = 25; observer_y = 0; observer_x = 50;
k = 1; % omega = k
num_rays = 100;
straight_path_time = sqrt((observer_x - source)^2 + (observer_y - source)^2) / k;
all_path=zeros(1,num_rays); all_delay = zeros(1,num_rays);

%for rays = 1:num_rays
%a = -pi/4;
%b = pi/4;

index=1;

for rays = linspace(-pi/4,pi/4,num_rays)
phi=rays;
r01 = screen_dist/cos(phi); r01_phase = exp(1i*k*r01);
r12 = sqrt((observer_x - screen_dist)^2 + (r01*sin(phi)-observer_y)^2); r12_phase = exp(1i*k*r12);
total_phase = r01_phase*r12_phase;

all_path(1,index)=total_phase;

time_of_path = (r01/k) + (r12/k);
time_delay = time_of_path - straight_path_time;
all_delay(1,index)=time_delay;

index=index+1;
end

secondary_paths = zeros(num_rays,num_rays); secondary_delays = zeros(num_rays,num_rays);
for i = 1:length(all_path)
% Replicate array at time delay to 
delay_at_step = all_delay(i);
path_at_step = all_path(i);
secondary_delays(i,1:num_rays)=repelem(delay_at_step,num_rays);
path_diff_at_step = repelem(path_at_step,num_rays) - all_path;
secondary_paths(i,1:num_rays)=path_diff_at_step;
end

secondary_paths = reshape(secondary_paths',[1,(num_rays^2)]);
secondary_delays = reshape(secondary_delays',[1,(num_rays^2)]);

scatter(secondary_delays,log(secondary_paths))
%scatter(log10(secondary_paths),log10(secondary_delays))

%scatter(all_path,all_delay)
%scatter(log(all_path),all_delay)
%scatter(log(all_delay),log(all_path))
% index=1;
% for rays = linspace(-pi/4,pi/4,num_rays)
% %phi = (b-a).*rand(1,1) + a;
% phi=rays;
% r01 = screen_dist/cos(phi); %r01_phase = exp(1i*k*r01);
% r12 = sqrt((observer_x - screen_dist)^2 + (r01*sin(phi)-observer_y)^2); %r12_phase = exp(1i*k*r12);
% total_path = r01 + r12;
% all_path(1,index)=total_path;
% 
% time_of_path = (r01/k) + (r12/k);
% time_delay = time_of_path - straight_path_time;
% all_delay(1,index)=time_delay;
% 
% index=index+1;
% end
% 
% all_phase=zeros(1,num_rays);
% for paths = 1:length(all_path)
%     all_phase(1,paths)=exp(1i*k*((2*all_path(paths))-sum(all_path)));
% end
% 
% 
% scatter(all_phase,all_delay)
