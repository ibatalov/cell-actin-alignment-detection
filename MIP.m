% Creates maximum intensity projection (MIP) from a z-stack
% input:
%       z_stack - z-stack to do the MIP from
%       numberOfChannels - total number of channels in the z-stack
%       channelsToAnalyze - a vector containing numbers of channels that
%       need to be present in the ouput

% output: a cell array containing the MIPs for each requested channel in
% the order they were requested
function [ max_projections ] = MIP(z_stack, numberOfChannels, channelsToAnalyze)
   
    numberOfSlices = size(z_stack, 1)/numberOfChannels;
    numberOfOutputChannels = length(channelsToAnalyze);
    max_projections = cell(numberOfOutputChannels, 1);
    
    for output_channel = 1 : numberOfOutputChannels
        channel = channelsToAnalyze(output_channel);
        max_projections{output_channel} = zeros(size(z_stack{channel}));
    end
    
    for slice = 1 : numberOfSlices
        for output_channel = 1 : numberOfOutputChannels
            channel = channelsToAnalyze(output_channel);
            max_projections{output_channel} = max(max_projections{output_channel}, double(z_stack{channel + (slice - 1)*numberOfChannels}));
        end
    end
end

