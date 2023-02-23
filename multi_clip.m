clear; clc; close;
directory_info = get_directory_info();
addpath(genpath(directory_info.chronux_folder));

%% Set wheither to group genotypes or not
% This script has two functions: a first pass look at all of the
% recordings individually, and also to compare the wildtype (wt) and mutant
% (mt) recordings. The `compare_genotypes` variable switches between these
% roles. If `compare_genotypes == false`, the plots ignore genotype data.
% If it's set to `true`, they will usually produce similar images, but they
% will annotate wich components have wich genotype.

compare_genotypes = true;


wt_color = 	"#0072BD";
mt_color = 	"#D95319";

%% Get clip information

all_clips = get_clip_metadata(); 
% Here clips mean distinct segments of recorded time; each clip has fields
% for many of its experimental conditions (like animal, genotype, and 
% temperature) so they can be sorted (or indexed to produce `clips` below). 
% They also have the location of the file to extract the clip out of. Note 
% that there can be multiple clips from one file.

clips = all_clips(all_clips.Temp == "37" & all_clips.Seizure == 0,:); 
% clips = all_clips(all_clips.Seizure == 1,:); 
% clips = all_clips(all_clips.Animal == "O" | all_clips.Seizure == 1,:); 


% `clips` contains the subset of the clip information we're interested in;
% for example, only the 37 degree data. Be careful modifying `clips` after 
% it's been created, because it's assumed the columns of `A` match up to 
% the rows of `clips`.


%% Construct A

% `A` is the data matrix extracted from the clips. Each column corresponds
% to one clip, and each row is a timepoint. So `A(:,1)` would be the the
% data extracted from the first clip.

% Fs stands for sampling frequency
% TODO: automatically check that Fs is accurate for each file
Fs = 4096;
downsampling_factor = 4;


% This determines the longest common sample length (for the height
% of A) and then pre-allocates A.
common_length = min(clips.Range(:,2) - clips.Range(:,1))*Fs;
A = nan([round(common_length/downsampling_factor), size(clips,1)]);
off_A = nan([round(common_length/downsampling_factor), size(clips,1)]);
% `off_A` is the matrix of off-channels; it can be dropped in anywhere `A`
% is; that might make it useful for comparisons?
% TODO: try using off_A for all these plots



for idx = 1:size(clips,1)
    c = get_lfp(clips.Filename(idx));
    channel = clips.("Better Channel")(idx);
    range = (clips.Range(idx,1) * Fs + 1);
    range = round(range):round(range+common_length-1);

    % Jansen21 downsamples as well; it's good for anti-ailiasing
    A(:,idx) = decimate(c(range, channel),downsampling_factor);
    A(:,idx) = detrend(A(:,idx));
    off_A(:,idx) = decimate(c(range, 3-channel),downsampling_factor);
    off_A(:,idx) = detrend(off_A(:,idx));
end

Fs = Fs/downsampling_factor;
time = (1:size(A,1))/Fs;


%% Plot all raw traces on top of each other
% this is more of a debugging cell; I would use the stack below for
% visualizaiton

figure;
if compare_genotypes
    plot(time, A(:,clips.Wildtype == 1), 'Color', wt_color); hold on;
    plot(time,A(:,clips.Wildtype == 0), 'Color', mt_color); hold off;
else
    plot(time, A(:,:));
end
xlabel("time (s)")
ylabel("EEG amplitude (µV)")
title("Downsampled EEG traces");

legend(clips.DisplayName)


%% Plot all raw traces (scaled) in a stack
% TODO: the traces look different in this plot; maybe explore that?

figure;
hold on;
for idx=1:size(A,2)
    if compare_genotypes
        if clips.Wildtype(idx) == 1
            plot(time, rescale(A(:,idx)) + idx-1, 'Color', wt_color);
        else
            plot(time, rescale(A(:,idx)) + idx-1, 'Color',mt_color); 
        end
    else
        plot(time, rescale(A(:,idx)) + idx-1); % note that the data is rescaled
    end
    
end
hold off;

if compare_genotypes
    legend('');
    hold on;
    plot([NaN NaN], [NaN NaN], 'Color', wt_color, 'DisplayName', 'Wildtype');
    plot([NaN NaN], [NaN NaN], 'Color', mt_color, 'DisplayName', 'Mutant');
    hold off;
    legend()
else
    legend(clips.DisplayName);
end

xlabel("time (s)")
ylabel("EEG amplitude (scaled)")
title("Downsampled EEG traces");

%% Show PSD Welch estimation
% This section estimates the power spectra of the signals; this essentially
% characterizes how much different frequencies of oscillation contribute to
% the signal.

% TODO: some other papers seem to get steeper curves than we seem to. Why?

window_length = round(1*Fs);
overlap_ratio = .5;
overlap = round(window_length * overlap_ratio);


[Pxx, f] = pwelch(A, window_length, overlap, [], Fs);
if compare_genotypes
    for i=1:size(Pxx,2)
        if compare_genotypes
            color = mt_color;
            if clips.Wildtype(i)
                color = wt_color;
            end
        end
        plot(f,log(Pxx(:,i)), 'Color', color); hold on;
    end
else
    plot(f,log(Pxx)); 
end
hold off;

xlabel("Frequency (Hz)")
ylabel("Power/frequency (dB/Hz)")
xlim([0,200])
legend(clips.DisplayName)



%% Plot some 60 Hz harmonics
% This is probably a dead-end analysis, but it's something we looked at, so
% I'll keep it. In the preliminary analysis, we thought that the odd(even?)
% harmonics of the 60Hz power-line noise had bigger peaks in the mutant
% animals. To explore this more, this cell plots the harmonics of 60Hz 
% noise for each clip.
% These peaks are quantified in the next cell, and the coloring of the
% plots produced here gives the key to that quantificaiton. I measured peak
% heights as standard deviations from baselinel and the baseline sample is 
% in yellow, the peak is plotted in orange, and ignored values are in blue.

mode = "mt";
% switching this between "welch" and something else toggles the way the
% power spectrum is calculated

if strcmp(mode, "welch")
    [Pxx, f] = pwelch(A, window_length, overlap, [], Fs);
    pxx = log(Pxx);

    center_r = 2;
    surround_r = 20;
    donut_hole_r = center_r;
else
    clear("params");
    params.Fs = Fs;
    NW = 10; % time-bandwith product
    params.tapers = [NW, round(2*NW-1)];
    
    [S,f] = mtspectrumc(A, params);
    pxx = log(S);

    
    center_r = 7;
    surround_r = 150;
    donut_hole_r = center_r + 50;
end


% here you can choose which harmonics to plot (e.g. only even or odd)
harmonics = 2:1:8;
heights = nan(numel(harmonics),size(A,2)); % this is the statistic plotted next cell
should_plot = true;

for i = 1:numel(harmonics)
    frequency = 60 * harmonics(i);
    [~,index] = min(abs(f - frequency));



    surround_slice = -surround_r:surround_r;
    
    center_slice = -center_r:center_r;
    

    donut_slice = [-surround_r:-donut_hole_r, donut_hole_r:surround_r];
    baseline = mean(pxx(donut_slice + index, :));
    baseline_std = std(pxx(donut_slice + index, :));
   
    heights(i,:) = mean((pxx(center_slice + index, :) - baseline)./baseline_std,1);

    if should_plot
        figure;

        for j = 1:size(A,2)
            subplot(1,size(A,2), j);
            plot(f((surround_slice + index)),pxx(surround_slice + index, j)); hold on;
            plot(f((center_slice + index)),pxx(center_slice + index, j));

            colors = get(gca,'colororder');
            plot(f((donut_slice(1:end/2) + index)),pxx(donut_slice(1:end/2) + index, j), "Color", colors(3,:));
            plot(f((donut_slice(end/2 +1:end) + index)),pxx(donut_slice(end/2 + 1:end) + index, j), 'c', "Color", colors(3,:));


            yline(baseline(j));
            yline(baseline(j) + baseline_std(j), '--');
            yline(baseline(j) - baseline_std(j), '--');
            if compare_genotypes
                star = "";
                if ~(clips.Wildtype(j) == 1)
                    star = "*";
                end
                title(clips.Animal(j) + star);

            else
                title(clips.Animal(j));
            end
            
        end

    end

end

%% Compare 60Hz harmonic statistics
% This cell quantifies the results of the preivous cell; it's pretty much
% only useful with `compare_genotypes` set to `true`.

figure;
heights(heights <0) = 0;
if compare_genotypes
    plot(harmonics,heights(:, clips.Wildtype==1),'.', 'color',wt_color, 'MarkerSize', 15); hold on;
    plot(harmonics,heights(:, clips.Wildtype==0),'.', 'color',mt_color, 'MarkerSize', 15); hold off;
else
    plot(harmonics,heights, '.', 'MarkerSize',15);

end
legend(clips.DisplayName)

%% Calculate comodulograms for all channels
% this calculates comodulograms for all of the samples in `A`/`clips`

iqr_normalization = false; 


% I sometimes use "carrier" for the frequency you take the phase of
% (usually the lower) and "modulated" for the frequency you take the
% amplitude of (usually the higher).
carrier_frequency_bins = linspace(2,14,1+12*2);
modulated_frequency_bins = linspace(10,200,1+19*2);
time_to_calculate_on = 10; % in seconds
% I calculate on a limited span of time because it speeds up the
% calculation

%TODO: check if there's a better way to set the time bins
n_time_bins = round(round(time_to_calculate_on*carrier_frequency_bins(1))*.75);
comodulograms = nan(numel(modulated_frequency_bins)-1, numel(carrier_frequency_bins)-1,size(A,2));

rcfb = [carrier_frequency_bins(1:end-1);carrier_frequency_bins(2:end)]';
rmfb = [modulated_frequency_bins(1:end-1);modulated_frequency_bins(2:end)]';
for k=1:size(comodulograms,3)
    if iqr_normalization
        [comodulograms(:,:,k), ~] = calculate_comodulogram2(A(1:time_to_calculate_on*Fs,k), rcfb, rmfb, Fs, 18,n_time_bins);
    else
        comodulograms(:,:,k) = calculate_comodulogram(A(1:time_to_calculate_on*Fs,k), rcfb, rmfb, Fs, 18);
    end

    
    fprintf("Done with comodulogram %d of %d.\n", k, size(comodulograms,3));
end

% saving these variables means you can load a comodulogram even if the
% `clips` variable has changed
comodulogram_genotypes = clips.Wildtype;
comodulogram_display_names = clips.DisplayName;

% calculating the comodulogram can take a long time, so this saves them to
% be loaded again if needed
filename = [directory_info.output_folder filesep "%s comodulogram.mat"];
filename = sprintf(filename,datestr(now,'yyyy-mm-dd_HH-MM'));
save(filename,'comodulograms', 'carrier_frequency_bins', 'modulated_frequency_bins', 'time_to_calculate_on', 'comodulogram_genotypes', 'comodulogram_display_names', 'iqr_normalization')

%% show comodulograms
% filename = [directory_info.output_folder filesep "2023-02-13_15-35 comodulogram.mat"];


% assuming this is run right after the last block, it loads the last
% generated file.
load(filename)


mod_freq_bin_centers = conv(modulated_frequency_bins, [.5, .5], 'valid');
carr_freq_bin_centers = conv(carrier_frequency_bins, [.5, .5], 'valid');


transformed_comodulograms = log(comodulograms);
crange = [min(transformed_comodulograms(:)), max(transformed_comodulograms(:))];
% crange = [0 10] * 10^-4;
for k=1:size(comodulograms,3)
    figure;
    imagesc(carr_freq_bin_centers,mod_freq_bin_centers, transformed_comodulograms(:,:,k), crange); 
    axis xy; 
    colorbar;
    added_text = "";
    if compare_genotypes
        if comodulogram_genotypes(k) == 1
            added_text = " (WT)";
        else
            added_text = " (MT)";
        end
    end
    title(comodulogram_display_names(k) + added_text)
    xlabel("frequency band for phase")
    ylabel("frequency band for amplitude")
end

%% Compare average comodulogram across WT and MT

wt_comod = mean(transformed_comodulograms(:,:,clips.Wildtype==1),3);
mt_comod = mean(transformed_comodulograms(:,:,clips.Wildtype==0),3);
genotype_difference = mt_comod - wt_comod;

figure;
imagesc(carr_freq_bin_centers,mod_freq_bin_centers, genotype_difference); 
axis xy; 
colorbar;
title("mean MT - WT comodulogram");
xlabel("frequency band for phase")
ylabel("frequency band for amplitude")

%% Calculate spectrograms for all channels

clear("params");
params.Fs = Fs; % units are Hz
% params.fpass = [0, params.Fs/2];
% params.fpass = [0, 200];
% params.pad = -1;
% params.err = [0, 0.05]; % 2 is jacknife

NW = 3;
params.tapers = [NW, round(2*NW-1)];


window_size = 1;
window_ratio = .1;
movingwin = [window_size, window_size*window_ratio];

[S,t,f] = mtspecgramc(A, movingwin, params);
%% Show spectrograms for all channels

for idx = 1:size(A,2)
    figure;
    plot_matrix(S(:,:,idx),t,f);
    title(clips.DisplayName(idx));
    ylim([0,200])
end

%% calculate coherency
% TODO: there might be someting here; two of the WT segretate at high
% frequency values
clear("params");
params.Fs = Fs;
% params.fpass = [0, 200];
NW = 20; % time-bandwith product
params.tapers = [NW, round(2*NW-1)];


[C, phi, S12,S1,S2, f] = coherencyc(A, off_A, params);
figure;
C = log(C);

if compare_genotypes
    plot(f,C(:,clips.Wildtype==1), 'Color', wt_color); hold on;
    plot(f,C(:,clips.Wildtype==0), 'Color', mt_color); hold off;
else
    plot(f,C);
end

xlim([0,200])

% figure;
% plot_vector(S1,f);



