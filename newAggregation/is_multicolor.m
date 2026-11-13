function tf = is_multicolor(SI_Info)
%IS_MULTICOLOR True if the TIFFs hold more than one imaging channel.
%
%   tf = is_multicolor(SI_Info)
%
%   Decided from the channels SAVED to the TIFFs (SI_Info.hChannels.channelSave), because that is what the
%   suite2p data contains. hChannels.channelsActive also lists channels that are only displayed: in
%   AnimalRB23_261022_1908 channelsActive = [1 2] but channelSave = 1 (one channel in the TIFFs, suite2p
%   nchannels = 1), which used to be taken for a two-colour recording and failed with
%   'Unrecognized field name "meanImg_chan2"'. channelsActive is only used when channelSave is missing.

tf = false;
if isempty(SI_Info) || ~isstruct(SI_Info) || ~isfield(SI_Info, 'hChannels')
    return
end
ch = SI_Info.hChannels;
if isfield(ch, 'channelSave') && ~isempty(ch.channelSave)
    tf = numel(ch.channelSave) > 1;
elseif isfield(ch, 'channelsActive')
    tf = numel(ch.channelsActive) > 1;
end
end
