% Example code demonstrating the calculation of the free-field and tank
% received time series at a receiver expected based on the generation of
% several different source wavefourms at a source.

% Written by Hayden Johnson, 2024-03-11

% specify constants
c = 1490; % speed of sound in m/s
beta_wall = -0.9; % wall and bottom reflection coefficient
beta_surface = -0.9; % surface reflection coefficient
cutoff_time = 10e-3; % minimum time for which all reflected paths are included, in s

% define tank size
Lx = 0.57;
Ly = 0.34;
Lz = 0.4;

% specify source position
x_source = 0.19;
y_source = 0.17;
z_source = 0.2;
r_source = [x_source; y_source; z_source];

% specify receiver_position
x_receiver = 0.36;
y_receiver = 0.17;
z_receiver = 0.11;
r_receiver = [x_receiver; y_receiver; z_receiver];

% create an array of several source waveforms
dt = 1/96000;
t = (0:dt:0.03-dt)';
f_source = 1e3*[1 2 4 8];
t0 = 0.005;
p_source = exp((t0-t)./(1e-3)).*sin(2*pi*f_source.*(t-t0));
p_source(t<t0,:) = 0;

% compute free-field received signal
p_receiver_free = compute_time_series_free_field(p_source,dt,r_source,r_receiver,c);

% compute tank received signal
p_receiver_tank = compute_time_series_with_tank_wall_reflection(p_source,dt,r_source,r_receiver,Lx,Ly,Lz,c,beta_wall,beta_surface,cutoff_time);

% plot time series
for i = 1:size(p_source,2)
    figure(i);
    clf;
    subplot(2,1,1);
    plot(t*1e3,real(p_source(:,i)));
    xlabel('Time (ms)');
    title('Source waveform');

    subplot(2,1,2);
    plot(t*1e3,real(p_receiver_free(:,i)),'displayname','Free field');
    hold on;
    plot(t*1e3,real(p_receiver_tank(:,i)),'displayname','Tank');
    hold off;
    xlabel('Time (ms)');
    title('Tank receiver waveforms');
    legend('location','northeast');
end