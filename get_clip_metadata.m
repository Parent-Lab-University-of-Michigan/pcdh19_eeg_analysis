function all_clips = get_clip_metadata()
% This file keeps track of the clips you've annotated from your EDF files.
% Its central object is the table `clips`, which annotates one clip per
% row. Note that there can be multiple clips from the same file (or even
% overlapping clips.) Usually when `clips` is used in a file, there are two
% copes: `all_clips` and `clips`. `all_clips` keeps track of all of the
% clips returned from this function (so all that are known), while `clips`
% is a filtered subset of interest.
% 
% inputs: none
% outputs:
%     clips: a (nx8) table with the following columns:
%         Animal: (string) the name of the animal (usually a letter)
%         Temp: (string) the experimental temperature the animal was kept at during
%             the recording
%         Wildtype: (1 or 0) wheither the animal is wildtype 
%         Seizure: (1 or 0) wheither the clip contains a seizure 
%         Range: (1x2 numeric) the start and stop times of the clip of
%             interest.
%         Better Channel: which channel is better (1 for left 2 for right)
%         Display Name: the name used in figure legends
%         Filename: the file where this clip can be found

directory_info = get_directory_info();

all_clips = table();

% this section will change if you use a different structure or add files
%%%%%%%%%%% START

all_clips(end+1,:) = {"B", "37", 0, 0, [500 900]/4, 1};
all_clips(end+1,:) = {"C", "37", 1, 0, [0, 136]/4, 1}; 
all_clips(end+1,:) = {"D", "37", 0, 0, [35, 350]/4, 2};
all_clips(end+1,:) = {"E", "37", 1, 0, [0, 375]/4, 2};
all_clips(end+1,:) = {"F", "37", 1, 0, [.5, 34], 2};

all_clips(end+1,:) = {"L", "37", 0, 0, [1.4 92], 2};
all_clips(end+1,:) = {"L", "42", 0, 0, [1 65], 2};
all_clips(end+1,:) = {"N", "37", 0, 0, [11 347], 1};
all_clips(end+1,:) = {"N", "42", 0, 0, [19 88], 1};
all_clips(end+1,:) = {"O", "37", 0, 0, [25 89], 1};
all_clips(end+1,:) = {"T", "37", 0, 0, [5 260], 1};
all_clips(end+1,:) = {"T", "42", 0, 0, [11 107], 1};

all_clips(end+1,:) = {"O", "42", 0, 1, [1 340], 2};
all_clips(end+1,:) = {"F", "37", 0, 1, [31 56], 2};
% 42 actually means 42.5

%%%%%%%%%%% END

all_clips.Properties.VariableNames = {'Animal', 'Temp', 'Wildtype', 'Seizure', 'Range', 'Better Channel'};
all_clips.DisplayName = all_clips.Animal + " " + all_clips.Temp;
all_clips.DisplayName(all_clips.Seizure==1) = all_clips.DisplayName(all_clips.Seizure==1) + " sz";


nonseizure_prepath = [directory_info.data_folder filesep 'baseline' filesep];
all_clips.Filename(all_clips.Seizure==0) = nonseizure_prepath + all_clips.Animal(all_clips.Seizure==0) + "_"+ all_clips.Temp(all_clips.Seizure==0) +'.edf';

seizure_prepath = [directory_info.data_folder filesep 'seizure' filesep];
all_clips.Filename(all_clips.Seizure==1,:) = seizure_prepath + all_clips.Animal(all_clips.Seizure==1,:) + "_"+ all_clips.Temp(all_clips.Seizure==1,:) +'.edf';

% all_clips = sortrows(all_clips, 'Wildtype');
end

