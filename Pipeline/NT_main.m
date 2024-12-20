% This is the main function of the algorithm
% Input: kilosort cluster label, channel map, mean waveforms, 
% Output: Unit match assignment
% For more comparisons, users need to write their own loops
function NT_main(input, mwf1, mwf2) 

chan_pos = input.chan_pos;
chan_map = input.chan_map;
% if the number of channels in mw differ from chan_map, assume mw is the
% original data including all channels, and select only those included in
% the sort.
[nChanPos, ~] = size(chan_pos);
[~, nChanMW, ~] = size(mwf1);

if nChanPos < nChanMW
    mwf1 = mwf1(:,chan_map+1,:);
    mwf2 = mwf2(:,chan_map+1,:);
end

% Estimate location 
wf_metrics1 = wave_metrics(mwf1, input, input.data_path1); %col 9,10,11 = x,z,y
wf_metrics2 = wave_metrics(mwf2, input, input.data_path2);

% Estimate drift
output.threshold = input.threshold;
output.z_mode = 0;
output = create_EMD_input(input, output, wf_metrics1, wf_metrics2, mwf1, mwf2, 'pre'); 
output = EMD_unit_match(input,output, 'pre');
[output.diffZ,edges] = z_estimate(input);
output.z_mode = kernelModeEstimate(output.diffZ);

h = findobj('type', 'figure');
curr_fig_count = numel(h);
figure(curr_fig_count+1)
histogram(output.diffZ,edges);
title(sprintf('drift estimate = %.1f um', output.z_mode));

% EMD unit matching
output = create_EMD_input(input, output, wf_metrics1, wf_metrics2, mwf1, mwf2, 'post'); 
output = EMD_unit_match(input, output, 'post');

% thresholding
output.results_wth = output.all_results_post(output.all_results_post(:,7) <= input.threshold,:); 
fprintf("    %d matches found \n", length(output.results_wth))

% Save
if ~exist(input.result_path, 'dir')
    mkdir(input.result_path);
end
save(fullfile(input.result_path, 'Input.mat'), "input");
save(fullfile(input.result_path, 'Output.mat'), "output");
end