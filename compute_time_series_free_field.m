function p_receiver = compute_time_series_free_field(p_source,dt,r_source,r_receiver,c)
    % Compute the time series observed at a receiver as a result of a time
    % series generated at a source, with free-field propagation. This 
    % function is essentially a wrapper for 
    % compute_free_field_greens_function.m to allow the user to work in the 
    % time domain instead of the frequency domain.
    
    % Inputs:
    % p_source: [NxM] Matrix of M length-N vectors specifying time series
    % at the source. Each column represents a source time series.
    % dt: Time step of the time series p_source (s)
    % r1: Vector position of the source (m) [3x1]
    % r2: Vector position of the receiver (m) [3x1]
    % c: Sound speed (m/s)
    
    % Outputs:
    % p_receiver: [NxM] Matrix of M length-N vectors which contain the
    % expected waveforms at the receiver resulting from the the time series
    % of p_source at the source.
    
    % Dependencies:
    % compute_free_field_greens_function.m
    
    % Written by Hayden Johnson, 2024-03-11

    %----------------------------------------------------------------------

    % compute frequencies of fft
    n = size(p_source,1);
    fft_freq = (0:n-1)'./(n*dt);
    fft_omega = 2*pi*fft_freq;
    
    % compute input signal fft
    fft_source = fft(p_source);
    
    % compute free-field greens function
    g_free = compute_free_field_greens_function(r_source,r_receiver,fft_omega,c);
    
    % compute free-field received signal
    fft_receiver = g_free.*fft_source;
    p_receiver = ifft(fft_receiver);
end