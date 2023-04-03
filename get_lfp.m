function [channels] = get_lfp(file, channel_prefix)
% returns the LFP signal stored in the given EDF file
%
% inputs: 
%     file: a filename(/path) specifying where to find the file of intrest
%     channel_prefix: a string (like "x1") giving the prefix of the channel
%         names in the edf file.
%
% outputs:
%     channels: an (nx2) array of LFP values, where the first column is the
%         left channel.

if nargin < 2
    channel_prefix = "x1";
end

data = edfread(file);
% TODO: this is orinally returned as a timetable; it might make more sense
% to keep them that way for long recordings? I'm not sure which functions
% would work on a timetable vs. matrix.

% % this gets metadata for the file
% edfinfo(file)


if ~any(strcmp(channel_prefix+"REF", data.Properties.VariableNames))
    error("Expected channel name not found in edf file. (see comment in code)");
    % if you're getting this error, you probably need to change the
    % channel_prefix variable. In my experience it's something like "x1" or
    % "x2", but you can check it by manually inspecting the EDF file. If
    % you're uncertain, run `data = edfread(file)` and then inspect the
    % `data` varaible in your workspace; look for a pattern like "x1REF",
    % and the "x1" will be your prefix.
end

reference_channel = cat(1,data.(channel_prefix+"REF"){:});



left_channel = cat(1,data.(channel_prefix+"LP"){:}) - reference_channel;
right_channel = cat(1,data.(channel_prefix+"RP"){:}) - reference_channel;
channels = [left_channel, right_channel];
end