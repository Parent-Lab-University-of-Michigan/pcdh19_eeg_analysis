clear; clc;
directory_info = get_directory_info();
addpath(genpath(directory_info.chronux_folder));


%% Creates a table of file information
all_clips = get_clip_metadata(); 
clips = all_clips(all_clips.Animal == "O",:); 

Fs = 4096;

%% compares left and right channels for each file
% this is where I decided the bounds to use in the clip metadata


for idx = 1:size(clips,1)
    c = get_lfp(clips.Filename(idx));

    if clips.("Better Channel")(idx) == 1
        l_color = "#D95319";
        r_color = "#EDB120";
    else
        l_color = "#EDB120";
        r_color = "#D95319";
    end
    
    axes = [];
    figure;

    subplot(2,1,1)
    time = (1:size(c,1))/Fs;
    plot(time,c(:,1)); hold on;

    l = clips.Range(idx,1);
    r = clips.Range(idx,2);

    [~, l] = min(abs(time-l));
    [~, r] = min(abs(time-r));

    plot(time(l:r), c(l:r,1), 'Color',l_color); hold off;

    xlabel('time (s)');
    ylabel('voltage (uV)');
    title(clips.DisplayName(idx) + " left channel")
    axis(gca, 'tight');
    axes(1) = gca;

    subplot(2,1,2)
    plot(time,c(:,2)); hold on;

    l = clips.Range(idx,1);
    r = clips.Range(idx,2);

    [~, l] = min(abs(time-l));
    [~, r] = min(abs(time-r));


    plot(time(l:r), c(l:r,2), 'Color',r_color); hold off;
    xlabel('time (s)');
    ylabel('voltage (uV)');
    title(clips.DisplayName(idx) + " right channel")
    axis(gca, 'tight');
    axes(2) = gca;
    

    linkaxes(axes, 'xy');
%     xlim([10,60])
%     ylim([-500, 500]);
end

%% This selects one clip in particular to look at for the rest of the file

clip_idx = 1;
c = get_lfp(clips.Filename(clip_idx));
time = (1:size(c,1))/Fs;


time1 = clips.Range(clip_idx,1);
time2 = clips.Range(clip_idx,2);

waveform = detrend(c((time1*Fs + 1):time2*Fs,1));
waveform_time = time((time1*Fs + 1):time2*Fs);

% c is the clip; it has a left and right channel

%% plots both channels as traces on the same plot

figure(101);
plot(time,c)
legend("left", "right")
ylabel("eeg (uV)")
xlabel("time (s)")
title("Comparison of left and right channels")

%% Welch periodogram

% This calculates the Welch periodogram for the marked time window for the
% current clip; this is usually selected to be the interval with the least
% noise. You can find the clip limits in `get_clip_metadata.m`.

figure(102);

window_length = round(1*Fs);
overlap_ratio = .5;
overlap = round(window_length * overlap_ratio);

[Pxx, f] = pwelch(waveform, window_length, overlap, [], 4096); 
plot(f, log(Pxx))
hold on;

xlabel("Frequency (Hz)")
ylabel("Power/frequency (dB/Hz)")
xlim([0,200])
legend(clips.DisplayName)

%% matlab Periodogram
% this is the same as the previous cell, but with the MATLAB builtin
% function

periodogram(waveform, [], [], Fs);

%% chronux periodogram
% this does the same thing as the previous two cells, but this cell uses
% the chronux library

clear("params");
params.Fs = Fs;
params.err = [2, 0.05]; % 2 is jacknife
% params.trialave = 1;

NW = 3; % TODO: this should be set systematically!
params.tapers = [ NW, round(2*NW-1)];

[S, f, Serr] = mtspectrumc(waveform, params);

figure(103);
plot_vector(S,f,[]); hold on;


%% walk through Harmonics of 60 Hz
% this is a precursor to the cell in the mutiple channel file that looks at
% the 60 hz harmonics; you can pretty much skip it

pxx = log(S);

n_harmonics = 5;
center_r = 75;
surround_r = 150;
heights = nan(n_harmonics,1);
should_plot = true;

for i = 1:n_harmonics
    frequency = 60 * i;
    [~,index] = min(abs(f - frequency));



    surround_slice = -surround_r:surround_r;
    
    center_slice = -center_r:center_r;
    

    donut_slice = [-surround_r:-center_r, center_r:surround_r];
    baseline = mean(pxx(donut_slice + index));
    baseline_std = std(pxx(donut_slice + index));

   
    heights(i) = sum((pxx(center_slice + index) - baseline)/baseline_std);

    if should_plot
        figure;
        plot(surround_slice,pxx(surround_slice + index)); hold on;
        plot(center_slice,pxx(center_slice + index));
        yline(baseline);
        yline(baseline + baseline_std, '--');
        yline(baseline - baseline_std, '--');

    end

end

% this is like the follow-up cell that plots the heights for each harmonic
figure(105);
plot(heights); hold on;


%% calculate taper spectrogram 
% this calculates a moving-time periodogram displayed in the cell below
% (using chronux)

clear("params");
params.Fs = Fs; % units are Hz
params.fpass = [0, 200];
NW = 3;
params.tapers = [ NW, 2*NW-1];

window_size = 1;
window_ratio = .5;
movingwin = [window_size, window_size*window_ratio];

[S,t,f] = mtspecgramc(waveform, movingwin, params);

%% show taper spectrogram
% this displays the spectrogram generated above

figure(104)
plot_matrix(S,t,f);

%% Shows spectral estimate with error bars

clear("params");
params.Fs = Fs; 
NW = 5; % TODO: this should be set systematically!
params.tapers = [ NW, round(2*NW-1)];

params.err = [2, 0.05]; % 2 is jacknife

[S, f, Serr] = mtspectrumc(waveform, params);

plot(f,log(S));
x = [f, fliplr(f)];
y = [log(Serr(2,:)), log(fliplr(Serr(1,:) ))];
p = patch(x, y, 'b');
p.FaceAlpha = .25;
p.EdgeAlpha = 0;

xlabel("Frequency")
ylabel("Power")

% plot_vector(S,f,[],Serr)