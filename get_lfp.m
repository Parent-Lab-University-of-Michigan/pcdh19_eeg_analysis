function [channels] = get_lfp(file)
% returns the LFP signal stored in the given EDF file
%
% inputs: 
%     file: a filename(/path) specifying where to find the file of intrest
%
% outputs:
%     channels: an (nx2) array of LFP values, where the first column is the
%         left channel.


data = edfread(file);
% TODO: this is orinally returned as a timetable; it might make more sense
% to keep them that way for long recordings? I'm not sure which functions
% would work on a timetable vs. matrix.

% % this gets metadata for the file
% edfinfo(file)

prefix = "x1";
if ~any(strcmp('x1RP', data.Properties.VariableNames))
    prefix="x2";
end
reference_channel = cat(1,data.(prefix+"REF"){:});



left_channel = cat(1,data.(prefix+"LP"){:}) - reference_channel;
right_channel = cat(1,data.(prefix+"RP"){:}) - reference_channel;
channels = [left_channel, right_channel];
end