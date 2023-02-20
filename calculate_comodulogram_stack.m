function [comodulograms, normalized_comodulogram] = calculate_comodulogram_stack(waveform, carrier_frequencies, modulated_frequencies, Fs, nbins, number_of_slices)
% This is similar to `calculate_comodulogram`, but it calculates a series
% of comodulograms based on windowing the input waveform into a sequence of
% equal-sized sub-waveforms. The idea is to determine how noisy each
% element of the comodulogram is; the returned comodulogram is calculated
% by a rank version of each cell's coefficient of variation. However, this
% could be used to also look at how comodulation changes over time. Note,
% however, that non-stationarity in the signal (like seizures) would make
% the returned `comodulogram` invalid.

% inputs:
%     waveform: the time series data as a vector
%     carrier_frequencies: a (Nx2) array, where each row is a frequency
%         range to get the angle from
%     modulated_frequencies: a (Mx2) array, where each row is a frequency
%         range to check for amplitude modulation
%     carrier_frequencies: the time series data as a vector
%     Fs: the sampling frequency of `waveform`
%     number_of_slices: the number of windows to break the signal into; be 
%         careful about breaking the waveform into windows smaller than the
%         period of the lowest frequency in `carrier_frequencies`.
% outputs:
%     comodulograms: a (MxNxn_zstack) array of comodulograms for each
%         window
%     normalized_comodulogram: a (MxN) array of modulation index values, 
%         where each row corresponds to a row of `modulated_frequencies` 
%         and each column corresponds to a row of `carrier_frequencies`.


c_phases = nan(numel(waveform), size(carrier_frequencies,1));

m_envelopes = nan(numel(waveform), size(modulated_frequencies,1));

for i = 1:size(carrier_frequencies)
    b = bandpass(waveform, carrier_frequencies(i,:), Fs);
    hb = hilbert(b);
    c_phases(:,i) = angle(hb);
end

for i = 1:size(modulated_frequencies)
%     b = bandpass(waveform, modulated_frequencies(i,:), Fs, ImpulseResponse="iir",Steepness=0.95);
    b = bandpass(waveform, modulated_frequencies(i,:), Fs, "Steepness",.94);
    hb = hilbert(b);
    m_envelopes(:,i) = abs(hb);
end


angle_bins = linspace(-pi,pi, nbins+1);
comodulograms = nan(size(m_envelopes,2), size(c_phases,2), number_of_slices);
time_bins = round(linspace(0,numel(waveform), number_of_slices+1));
time=(1:numel(waveform))';

for l=1:size(comodulograms,3)
    for k=1:size(comodulograms,2)
        for j=1:size(comodulograms,1)
            phase = c_phases(:,k);
            envelope = m_envelopes(:,j);
    
    
            lg_amplitude = nan(size(angle_bins,2)-1,1);
            
            for i=1:(numel(lg_amplitude))
                phase_mask = phase > angle_bins(i);
                if i ~= numel(angle_bins)
                    phase_mask = phase_mask & (phase <= angle_bins(i+1));
                end

                time_mask = time > time_bins(l);
                if i ~= numel(time_bins)
                    time_mask = time_mask & (time <= time_bins(l+1));
                end

                lg_amplitude(i) = mean(envelope(phase_mask & time_mask));

            end
            
            P = lg_amplitude / sum(lg_amplitude);
            
            U = ones(nbins,1)/nbins;
            KL = P.*(log(P)-log(U));
            KL = sum(KL);
            
            
            
            comodulograms(j,k,l) = KL/log(nbins);
        end
    end
end

normalized_comodulogram = median(comodulograms,3) ./ iqr(comodulograms,3);
end

