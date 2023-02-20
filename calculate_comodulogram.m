function comodulogram = calculate_comodulogram(waveform, carrier_frequencies, modulated_frequencies, Fs, nbins)
% This is equivalent to calling `get_modulation` on the same waveform over 
% a grid of frequencies, but caches data to speed the calculation up.

% inputs:
%     waveform: the time series data as a vector
%     carrier_frequencies: a (Nx2) array, where each row is a frequency
%         range to get the angle from
%     modulated_frequencies: a (Mx2) array, where each row is a frequency
%         range to check for amplitude modulation
%     carrier_frequencies: the time series data as a vector
%     Fs: the sampling frequency of `waveform`
% outputs:
%     comodulogram: a (MxN) array of modulation indices, where each row
%         corresponds to a row of `modulated_frequencies` and each column
%         corresponds to a row of `carrier_frequencies`.

% TODO: should there be overlaps between the frequency bins?

c_phases = nan(numel(waveform), size(carrier_frequencies,1));

m_envelopes = nan(numel(waveform), size(modulated_frequencies,1));

for i = 1:size(carrier_frequencies)
    b = bandpass(waveform, carrier_frequencies(i,:), Fs);
    hb = hilbert(b);
    c_phases(:,i) = angle(hb);
end

for i = 1:size(modulated_frequencies)
    b = bandpass(waveform, modulated_frequencies(i,:), Fs);
    hb = hilbert(b);
    m_envelopes(:,i) = abs(hb);
end


angle_bins = linspace(-pi,pi, nbins+1);
comodulogram = nan(size(m_envelopes,2), size(c_phases,2));

for k=1:size(comodulogram,2)
    for j=1:size(comodulogram,1)
        phase = c_phases(:,k);
        envelope = m_envelopes(:,j);


        lg_amplitude = nan(size(angle_bins,2)-1,1);
        
        for i=1:(numel(lg_amplitude))
            mask = phase > angle_bins(i);
            if i ~= numel(angle_bins)
                mask = mask & (phase <= angle_bins(i+1));
            end
            lg_amplitude(i) = mean(envelope(mask));
        end
        
        P = lg_amplitude / sum(lg_amplitude);
        
        U = ones(nbins,1)/nbins;
        KL = P.*(log(P)-log(U));
        KL = sum(KL);
        
        
        
        comodulogram(j,k) = KL/log(nbins);
    end
end

end

