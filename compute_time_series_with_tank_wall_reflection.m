function p_receiver = compute_time_series_with_tank_wall_reflection(p_source,dt,r_source,r_receiver,Lx,Ly,Lz,c,beta_wall,beta_surface,cutoff_time)
    % Compute the time series observed at a receiver as a result of a time
    % series generated at a source in a tank. This function is essentially 
    % a wrapper for compute_tank_greens_function.m to allow the user to
    % work in the time domain instead of the frequency domain.
    
    % Inputs:
    % p_source: [NxM] Matrix of M length-N vectors specifying time series
    % at the source. Each column represents a source time series.
    % dt: Time step of the time series p_source (s)
    % r_source: Vector position of the source (m) [3x1]
    % r_receiver: Vector position of the receiver (m) [3x1]
    % Lx, Ly, Lz: Dimensions of the tank in each coordinate (m)
    % c: Sound speed (m/s)
    % beta_wall: Reflection coefficient for the 5 non-surface walls of the
    % tank
    % beta_surface: Reflection coefficient for the water surface
    % cutoff_time: Time over which to sum reflected paths (s) 
    
    % Outputs:
    % p_receiver: [NxM] Matrix of M length-N vectors which contain the
    % expected waveforms at the receiver resulting from the the time series
    % of p_source at the source.
    
    % Source and receiver positions are specified in meters in a coordinate
    % system in which the origin lies at one of the vertices of the tank,
    % and x, y, and z are all increasing into the tank, such that the tank
    % walls lie at x=0, x=Lx, y=0, y=Ly, z=0, and the surface is at z=Lz.
    
    % Every reflected path which arrives before the cutoff time is included
    % in the summation, with time t=0 corresponding to the start of the 
    % signal production at the source. Note that some paths which arrive 
    % after the cutoff time will also be included in the summation, so the
    % output is no longer exact after the cutoff time. Ideally, the cutoff
    % time should be chosen to be sufficiently long such that reflected
    % paths which take longer than this to arrive have been sufficiently
    % dampened as to not contribute meaningfully to the signal measured at
    % the receiver.
    
    % Dependencies:
    % compute_tank_greens_function.m
    % compute_free_field_greens_function.m
    
    % Written by Hayden Johnson, 2024-03-11
    
    %----------------------------------------------------------------------
    
    % compute frequencies of fft
    n = size(p_source,1);
    fft_freq = (0:n-1)'./(n*dt);
    fft_omega = 2*pi*fft_freq;
    
    % compute input signal fft
    fft_source = fft(p_source);
    
    % compute greens function
    g_tank = compute_tank_greens_function(r_source,r_receiver,fft_omega,Lx,Ly,Lz,c,beta_wall,beta_surface,cutoff_time);

    % compute received signal
    fft_receiver_tank = g_tank.*fft_source;
    p_receiver = ifft(fft_receiver_tank);
end