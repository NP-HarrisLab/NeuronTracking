% Main file to run to track units, with reference data
% Adapted from Yuan's Example_run.m
ks_folders = ["D:\T02\catgt_20240705_Bank0_SemiCheck_R_g0\20240705_Bank0_SemiCheck_R_g0_imec0\imec0_ks25", "D:\T02\catgt_20240710_T02_Bank0_SemiCheck_R_g0\20240710_T02_Bank0_SemiCheck_R_g0_imec0\imec0_ks25"];
npy_matlab_path = "C:\Users\tdeweese\Documents\SpikeSorting\npy-matlab";
output_path = "D:\tracking";
run_tracking(ks_folders, npy_matlab_path, output_path)
function run_tracking(ks_folders, npy_matlab_path, output_path, mean_wf_name, metrics_name)
    arguments
        ks_folders (1,:) string {mustBeLongerThanOne(ks_folders)}
        npy_matlab_path (1,:) string
        output_path (1,:) string
        mean_wf_name (1,:) string {mustBeFileType(mean_wf_name, ".npy")} = "mean_waveforms.npy"
        metrics_name (1,:) string {mustBeFileType(metrics_name, ".csv")} = "metrics.csv"
    end

    % Arguments
    %   ks_folders: list of kilosort folders to process in order.
    %        Each ks output must include:
    %               channel_map.npy
    %               channel_positions.npy
    %               mean_waveforms.npy (if other name, specify in
    %               mean_wf_name)
    %               metrics.csv (if other name, specify in metrics_name)
    %                   just needs 'cluster_id' and 'firing_rate' 
    %               cluster_group.tsv
    %   npy_matlab_path: path to the npy-matlab package
    %   mean_wf_name: mean waveforms filename. Default: mean_waveforms.npy
    %   metrics_name: metrics filename. Default: metrics.csv

    % add packages
    neuron_tracking_path = fileparts(mfilename('fullpath'));
    addpath(genpath(neuron_tracking_path))
    addpath(genpath(npy_matlab_path))
    
    %----------EDIT to define paths for input----------
    input.input_paths = ks_folders;
    input.output_path = output_path;
    input.EMD_path = fullfile(input.output_path,'EMD_input\'); % Directory for this pipeline to store input, create before running
    input.wf_name = mean_wf_name;
    input.KSLabel_name = 'cluster_group.tsv'; % name for files of unit calls
    input.chan_pos_name = 'channel_positions.npy';
    input.chan_map_name = 'channel_map.npy'; 
    input.shank = -1; % To include all shanks, set to -1, for 2.0 probes, set to 0-based probe index to limit to single shank
    
    
    % Input data characteristics
    input.fs = 30000; %acquisition rate, Neuropixels default
    input.ts = 82; %wf time samples, set by C_Waves
    
    % parameters for running match
    input.l2_weights = 1500;
    input.threshold = 10; %z distance threshold for matched units, 10um recommended, can change if needed
    input.validation = 0; % set to 1 if you are including validation data, see README for format
    input.xStep = 32; %space between columns of sites, um (NP 2.0 = 32, can find in the channel map)
    input.zStep = 15; %space between rows of sites, um (NP 2.0 = 15, can find in the channel map)
    input.dim_mask = logical([1,1,1,0,0,0,0,0,0,1]); %default = x,z,y position, waveform distance
    input.dim_mask_physical = logical([1,1,1,0,0,0,0,0,0,0]);
    input.dim_mask_wf = logical([0,0,0,0,0,0,0,0,0,1]);
    input.diagDistCalc = true; % set to true to include separate calcualtion of distance and waveform sim matrices
    
    numData = length(input.input_paths); %Number of datasets to match
    
    
    
    %----------Unit tracking----------
    if exist(input.EMD_path, 'dir') == 0
        mkdir(input.EMD_path);
    end
    
    % Read in chan_pos and chan_map from first day; all other days should match
    day1 = input.input_paths(1);
    chan_pos = readNPY(fullfile(day1, input.chan_pos_name));
    chan_map = readNPY(fullfile(day1, input.chan_map_name));
    
    % Find match of all datasets 
    for id = 1:numData-1
        input.data_path1 = input.input_paths(id); % first dataset, NEED CHANGE 
        input.data_path2 = input.input_paths(id+1); % second dataset, NEED CHANGE
        result_dir = sprintf('result_%d_%d',id,id+1);
        input.result_path = fullfile(input.output_path,result_dir); %result directory 
        input.input_name = ['input',num2str(id),'.mat']; 
        input.input_name_post = ['input_post',num2str(id),'.mat']; 
        input.filename_pre = ['EMD_pre',num2str(id),'.mat']; 
        input.filename_post = ['EMD_post',num2str(id),'.mat']; 
        input.chan_pos = chan_pos;
        input.chan_map = chan_map;
        mwf1 = readNPY(fullfile(input.data_path1, input.wf_name));
        mwf2 = readNPY(fullfile(input.data_path2, input.wf_name));
        NT_main(input, mwf1, mwf2);
    end
    
    
    
    %% ----------Plot matched units----------
    % Select a first datasets (integer in [1:numData-1] and
    plot_id = 1;
    % Select unit index in the list of matches with z distance < threshold
    plot_unit_index = 8;
    
    input_plot = input;
    input_plot.data_path1 = input.input_paths(plot_id);
    input_plot.data_path2 = input.input_paths(plot_id+1);
    result_dir = sprintf('result_%d_%d',plot_id,plot_id+1);
    input_plot.result_path = fullfile(input.output_path,result_dir); %result directory 
    input_plot.filename_post = ['EMD_post',num2str(plot_id),'.mat']; 
    
    % plot unit locations and matches of curr_id and curr_id+1
    % includes all matches (z distance threshold not applied).
    plot_unit(input_plot);
    
    % Plot sample waveform from 
    plot_out = load(fullfile(input_plot.result_path,'Output.mat'));
    % find matched units with zdist < 10 um
    all_matches = plot_out.output.all_results_post;
    thresh_match_ind = all_matches(:,7) < 10;
    thresh_matches = all_matches(thresh_match_ind,:);
    % in all_matches:
    %     column 2 = unit index for dataset (plot_id + 1)
    %     column 3 = unit index on dataset (plot_id)
    plot_waveform(input_plot, thresh_matches(plot_unit_index,3),thresh_matches(plot_unit_index,2));
    
    % Z distance distribution
    plot_z(fullfile(input_plot.result_path));
    
    
    % waveform vs physical distance. INPUT = result path 
    plot_dist(input_plot);
    
    
    %% Find chains
    all_input = cell(numData-1, 1);
    all_output = cell(numData-1, 1);
    for id = 1:numData-1 % Load data
        result_dir = sprintf('result_%d_%d',id,id+1);
        all_input{id} = load(fullfile(input.output_path,result_dir, "Input.mat")); 
        all_output{id} = load(fullfile(input.output_path,result_dir,'Output.mat'));
    end
    [chain_all,z_loc,len] = chain_summary(all_input,all_output,numData,input.output_path);
    fprintf('Completed chain search. \n')
    
    if ~isempty(len)
    
    %----------Plot chains of interest (waveform, firing rate, original location, drift-corrected location, L2)----------
    full_chain = chain_all(len == numData,:); %find chains with length across all datasets
    [L2_value,fr_all,fr_change,x_loc_all,z_loc_all] = chain_stats(all_input,all_output,full_chain,numData,input.input_path);
    
    numChain = size(full_chain,1);
    ichain = 1; %which chain to plot, please enter a number between 1 and numChain as input, NEED CHANGE  
    
    assert(ichain<numChain, "ichain must be less than %d", numChain)
    figure()
    for id = 1:numData-1
        % plot waveform
        plot_wf(all_input, full_chain, L2_value, chan_pos, numData, ichain, id);
    end
    % plot firing rate
    plot_fr(fr_all, fr_change, numData, ichain);
    % plot location
    plot_loc(all_input,x_loc_all,z_loc_all, chan_pos, numData,ichain)
    
    end
end

function mustBeFileType(filename, ext)
    if ~endsWith(filename, ext)
        error("The provided name ('%s') must end with ('%s')", filename, ext)
    end
end

function mustBeLongerThanOne(list)
    if length(list) <= 1
        error("Must have more than one ks_folder to perform tracking")
    end
end
