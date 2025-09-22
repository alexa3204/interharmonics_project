function V_pos_waveform = positive_sequence_waveform(Va, Vb, Vc, f, t)
% POSITIVE_SEQUENCE_WAVEFORM Calculate positive sequence waveform from three-phase waveforms
%
% Inputs:
%   Va, Vb, Vc - Three-phase waveforms as column vectors (same length)
%   f - Fundamental frequency in Hz
%   t - Time vector as column vector (same length as voltage vectors)
%
% Output:
%   V_pos_waveform - Positive sequence waveform as column vector
%
% Theory:
%   1. Extract phasor from each phase using correlation method
%   2. Apply symmetrical components transformation: V1 = (1/3)[Va + a*Vb + a²*Vc]
%   3. Reconstruct time-domain waveform from positive sequence phasor
%   where a = e^(j2π/3) is the 120° complex operator
%
% Example:
%   t = (0:0.0001:1/60)';  % One cycle at 60 Hz, column vector
%   Va = 100*cos(2*pi*60*t);
%   Vb = 100*cos(2*pi*60*t - 2*pi/3);
%   Vc = 100*cos(2*pi*60*t + 2*pi/3);
%   V_pos = positive_sequence_waveform(Va, Vb, Vc, 60, t);
%   plot(t, [Va, Vb, Vc, V_pos]);
%   legend('Va', 'Vb', 'Vc', 'V_{pos}');

    % Validate inputs
    if length(Va) ~= length(Vb) || length(Vb) ~= length(Vc) || length(Vc) ~= length(t)
        error('All input vectors must have the same length');
    end
    if f <= 0
        error('Frequency must be positive');
    end
    
    % Ensure column vectors
    Va = Va(:);
    Vb = Vb(:);
    Vc = Vc(:);
    t = t(:);
    
    % Define the complex operator 'a' for 120° phase shift
    a = exp(1j * 2 * pi / 3);  % a = e^(j2π/3)
    
    % Extract phasors from each phase using correlation method
    Va_phasor = extract_phasor(Va, f, t);
    Vb_phasor = extract_phasor(Vb, f, t);
    Vc_phasor = extract_phasor(Vc, f, t);
    
    % Calculate positive sequence phasor using symmetrical components
    V_pos_phasor = (1/3) * (Va_phasor + a * Vb_phasor + a^2 * Vc_phasor);
    
    % Reconstruct positive sequence waveform
    V_pos_waveform = real(V_pos_phasor * exp(1j * 2 * pi * f * t));
end

function phasor = extract_phasor(signal, f, t)
    % Extract phasor using correlation method (returns peak magnitude phasor)
    %
    % For a signal: signal = A*cos(2πft + φ)
    % The phasor representation is: A∠φ = A*e^(jφ)
    
    % Create reference signals
    cos_ref = cos(2*pi*f*t);
    sin_ref = sin(2*pi*f*t);
    
    % Calculate correlations
    % For signal = A*cos(2πft + φ) = A*cos(φ)*cos(2πft) - A*sin(φ)*sin(2πft)
    cos_corr = mean(signal .* cos_ref) * 2;  % = A*cos(φ)
    sin_corr = mean(signal .* sin_ref) * 2;  % = -A*sin(φ)
    
    % Extract magnitude and phase
    magnitude = sqrt(cos_corr^2 + sin_corr^2);
    phase_rad = atan2(-sin_corr, cos_corr);
    
    % Create phasor (peak magnitude)
    phasor = magnitude * exp(1j * phase_rad);
end