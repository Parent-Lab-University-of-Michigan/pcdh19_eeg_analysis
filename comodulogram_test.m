% This file walks through how to calculate modulation indexes on synthetic
% data. Near the end it also compares possible comodulogram plots.

clc; close all; clear;
directory_info = get_directory_info();

%% Create a synthetic signal, y, to analyse
Fs = 4096; % sampling frequency

T = 13.4; % duration of recording
N = T*Fs; % number of samples taken
t = linspace(0,T,N); % time for our signal


y_carrier = pi*sin(t * (2*pi) * 8); % the carrier signal
y_modulated = .5*sin(t * (2*pi) * 93) .* (1 + .45 * y_carrier); % the modulated signal
y_full = y_carrier+y_modulated;



y_full = y_full + .1*randn(size(y_full)); % random noise
% y = y + 1000; % DC offset

% this replicates our downsampling in the real data
downsampling_factor = 8;
y = decimate(y_full,downsampling_factor);

Fs = Fs/downsampling_factor;
old_t = t;
t = decimate(t,downsampling_factor);


%% plot y and its parts
figure;
tiledlayout(3,1);


ax1 = nexttile;
plot(old_t,y_carrier);
title("Carrier frequency");


ax2 = nexttile;
plot(old_t,y_modulated);
title("Modulated frequency");


ax3 = nexttile;
plot(old_t,y_full);
title("Composite signal");


linkaxes([ax1 ax2 ax3], 'x');
xlim([0 .5] + .5);

%% Compute bandpasses, envelope, and angles of y

carrier_bandpass = bandpass(y,[7, 9], Fs);
hilbert_carrier_bandpass = hilbert(carrier_bandpass);
phase = angle(hilbert_carrier_bandpass);

modulated_bandpass = bandpass(y,[85,95], Fs);
hilbert_modulated_bandpass = hilbert(modulated_bandpass);
amplitude = abs(hilbert_modulated_bandpass);
    

%% Work Backwards from downsampled y and compare to true values
figure;
tiledlayout(3,2);


ax1 = nexttile(1);
plot(old_t,y_carrier);
title("Carrier true signal");


ax2 = nexttile(3);
plot(old_t,y_modulated);
title("Modulated true signal");


ax3 = nexttile(5);
plot(old_t,y_full);
title("Composite signal");


ax4 = nexttile(6);
plot(t,y);
title("Downsampled Signal");


ax5 = nexttile(4);
plot(t,modulated_bandpass); hold on;
plot(t,amplitude);
title("Estimated modulated signal");



ax6 = nexttile(2);
plot(t,carrier_bandpass); hold on;
plot(t,phase);
title("Estimated carrier signal");

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x');

xlim([0 .5] + .5);

%% calculate mean modulated amplitudes for each phase

number_of_bins = 18;
angle_bins = linspace(-pi,pi,number_of_bins+1);
mean_amplitude_per_phase = nan(number_of_bins,1);

for idx=1:number_of_bins
    mask = phase > angle_bins(idx);
    if idx ~= number_of_bins
        mask = (mask) & (phase <= angle_bins(idx+1));
    end
    mean_amplitude_per_phase(idx) = mean((amplitude(mask)));
end


angle_centers = conv(angle_bins, [.5,.5], 'valid');
subplot(1,2,1);
polarplot(angle_centers, mean_amplitude_per_phase, '.', 'MarkerSize',40);

subplot(1,2,2)
bar(angle_centers, mean_amplitude_per_phase);

%% calculate modulation index

amplitude_distribution = mean_amplitude_per_phase / sum(mean_amplitude_per_phase);

uniform_distribution = ones(number_of_bins,1)/number_of_bins;
KL = amplitude_distribution.*(log(amplitude_distribution)-log(uniform_distribution));
KL = sum(KL);

modulation_index = KL/log(number_of_bins);
fprintf("Modulation is: %g\n", modulation_index);


%% calculate MI over a grid (a comodulogram)

carrier_frequency_bins = linspace(2,14,1+12*2);
modulated_frequency_bins = linspace(5,200,1+19.5*2);

rcfb = [carrier_frequency_bins(1:end-1);carrier_frequency_bins(2:end)]';
rmfb = [modulated_frequency_bins(1:end-1);modulated_frequency_bins(2:end)]';

slice_length = 10;
time_slice = 1:slice_length*Fs;
number_of_windows = round(round(slice_length*carrier_frequency_bins(1))*.75);

tic;
[comodulograms,~]= calculate_comodulogram_stack(y(time_slice), rcfb, rmfb, Fs, 18, number_of_windows);
toc;

%% Compare the normalized to unnormalized comodulograms

normalized_comodulogram = median(comodulograms,3) ./ iqr(comodulograms,3);

figure;
mod_freq_bin_centers = conv(modulated_frequency_bins, [.5, .5], 'valid');
carr_freq_bin_centers = conv(carrier_frequency_bins, [.5, .5], 'valid');

subplot(1,2,1)
imagesc(carr_freq_bin_centers,mod_freq_bin_centers, normalized_comodulogram); 
axis xy; 
colorbar;
title("Normalized comodulogram");

subplot(1,2,2)
imagesc(carr_freq_bin_centers,mod_freq_bin_centers, mean(comodulograms,3)); 
axis xy; 
colorbar;
title("Mean comodulogram");

%% This creates a video out of the comodulogram stack

figure;
im = imagesc(carr_freq_bin_centers,mod_freq_bin_centers, (comodulograms(:,:,1))); 
axis xy; 

writerObj = VideoWriter([directory_info.output_folder filesep 'comodulogram_stack.avi']); %// initialize the VideoWriter object
open(writerObj);
for t = 2:size(comodulograms,3)
   set(im,'CData',(comodulograms(:,:,t)));
   F = getframe;
   writeVideo(writerObj,F);
end
close(writerObj);


%% inspect the frequency domain of the signal (something in the process split up our peaks?)

window_length = round(1*Fs);
overlap_ratio = .5;
overlap = round(window_length * overlap_ratio);

pwelch(y, window_length, overlap, [], Fs);


