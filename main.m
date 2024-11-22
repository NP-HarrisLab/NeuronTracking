% Main file to run to track units, with reference data
% example with optional arguments:  main(ks_folders, npy_matlab_path, output_path, 'mean_wf_name', "ksproc_mean_waveforms.npy", 'cluster_group', "cluster_KSLabel.tsv");

% Arguments
%   ks_folders: list of kilosort folders to process in order.
%        Each ks output must include:
%               channel_map.npy
%               channel_positions.npy
%               mean_waveforms.npy (if other name, specify in
%               mean_wf_name)
%               metrics.csv (if other name, specify in metrics_name)
%                   needs columns 'cluster_id' and 'firing_rate'
%               cluster_group.tsv (if other name, specify in cluster_group)
%   npy_matlab_path: path to the npy-matlab package
%   mean_wf_name: mean waveforms filename. Default: mean_waveforms.npy
%   metrics_name: metrics filename. Default: metrics.csv
%   cluster_group: cluster group filename. Default: cluster_group.tsv
% arguments
%     ks_folders (1,:) string {mustBeLongerThanOne(ks_folders)}
%     npy_matlab_path (1,:) string
%     output_path (1,:) string
%     mean_wf_name (1,:) string {mustBeFileType(mean_wf_name, ".npy")} = "mean_waveforms.npy"
%     metrics_name (1,:) string {mustBeFileType(metrics_name, ".csv")} = "metrics.csv"
%     cluster_group(1,:) string {mustBeFileType(cluster_group, ".tsv")} = "cluster_group.tsv"
% end


function main(ks_folders, npy_matlab_path, output_path, varargin)
    %----------Parse inputs----------
    
    p = inputParser;
    addRequired(p, 'ks_folders', @(x) (isstring(x) || ischar(x)) && length(x) > 1);
    addRequired(p, 'npy_matlab_path', @(x) ischar(x) || isstring(x));
    addRequired(p, 'output_path', @(x) ischar(x) || isstring(x));
    
    % Optional arguments with default values
    addParameter(p, 'mean_wf_name', 'mean_waveforms.npy', @(x) mustBeFileType(x, ".npy"));
    addParameter(p, 'metrics_name', 'metrics.csv', @(x) mustBeFileType(x, ".csv"));
    addParameter(p, 'cluster_group', 'cluster_group.tsv', @(x) mustBeFileType(x, ".tsv"));
    
    parse(p, ks_folders, npy_matlab_path, output_path, varargin{:});
    
    ks_folders = p.Results.ks_folders;
    npy_matlab_path = p.Results.npy_matlab_path;
    output_path = p.Results.output_path;
    mean_wf_name = p.Results.mean_wf_name;
    metrics_name = p.Results.metrics_name;
    cluster_group = p.Results.cluster_group;
    
    % add packages
    neuron_tracking_path = fileparts(mfilename('fullpath'));
    addpath(genpath(neuron_tracking_path))
    addpath(genpath(npy_matlab_path))
    
    %----------EDIT to define paths for input----------
    input.input_paths = ks_folders;
    input.output_path = output_path;
    input.EMD_path = fullfile(input.output_path,'EMD_input\'); % Directory for this pipeline to store input, create before running
    input.wf_name = mean_wf_name;
    input.metrics_name = metrics_name;
    input.KSLabel_name = cluster_group; % name for files of unit calls
    input.chan_pos_name = 'channel_positions.npy';
    input.chan_map_name = 'channel_map.npy';
    input.shank = -1; % To include all shanks, set to -1, for 2.0 probes, set to 0-based probe index to limit to single shank
    
    
    % Input data characteristics
    input.fs = 30000; % default acquisition rate, will read in .meta if possible
    input.ts = 82; % wf time samples, set by C_Waves
    
    % parameters for running match
    input.l2_weights = 1500;
    input.threshold = 10; % z distance threshold for matched units, 10um recommended, can change if needed
    input.validation = 0; % set to 1 if you are including validation data, see README for format
    input.xStep = 32; % space between columns of sites, um (NP 2.0 = 32, can find in the channel map)
    input.zStep = 15; % space between rows of sites, um (NP 2.0 = 15, can find in the channel map)
    input.dim_mask = logical([1,1,1,0,0,0,0,0,0,1]); % default = x,z,y position, waveform distance
    input.dim_mask_physical = logical([1,1,1,0,0,0,0,0,0,0]);
    input.dim_mask_wf = logical([0,0,0,0,0,0,0,0,0,1]);
    input.diagDistCalc = true; % set to true to include separate calcualtion of distance and waveform sim matrices
    
    numData = length(input.input_paths); % Number of datasets to match
    
    
    
    %----------Unit tracking----------
    if exist(input.EMD_path, 'dir') == 0
        mkdir(input.EMD_path);
    end
    
    % Read in chan_pos and chan_map from first day; ensure all other days match
    day1 = input.input_paths(1);
    chan_pos = readNPY(fullfile(day1, input.chan_pos_name));
    chan_map = readNPY(fullfile(day1, input.chan_map_name));
    for id = 2:numData
        day = input.input_paths(id);
        chan_pos2 = readNPY(fullfile(day, input.chan_pos_name));
        chan_map2 = readNPY(fullfile(day, input.chan_map_name));
        assert(isequal(chan_pos, chan_pos2), "chan_pos must be the same across all datasets")
        assert(isequal(chan_map, chan_map2), "chan_map must be the same across all datasets")
    end
    
    % Find match of all datasets
    for id = 1:numData-1
        input.data_path1 = input.input_paths(id); % first dataset
        input.data_path2 = input.input_paths(id+1); % second dataset
        
        data_name1 = getDataName(input.data_path1);
        data_name2 = getDataName(input.data_path2);
        fprintf('Comparing datasets %s and %s \n', data_name1, data_name2)
        
        result_dir = sprintf('result_%d_%d', id,id+1);
        input.result_path = fullfile(input.output_path, result_dir); % result directory
        input.input_name = ['input', num2str(id), '.mat'];
        input.input_name_post = ['input_post',num2str(id), '.mat'];
        input.filename_pre = ['EMD_pre', num2str(id), '.mat'];
        input.filename_post = ['EMD_post', num2str(id), '.mat'];
        input.chan_pos = chan_pos;
        input.chan_map = chan_map;
        mwf1 = readNPY(fullfile(input.data_path1, input.wf_name));
        mwf2 = readNPY(fullfile(input.data_path2, input.wf_name));
        NT_main(input, mwf1, mwf2);
    end
    
    
    %----------Chain Summary----------
    all_input = cell(numData-1, 1);
    all_output = cell(numData-1, 1);
    for id = 1:numData-1 % Load data
        result_dir = sprintf('result_%d_%d',id,id+1);
        all_input{id} = load(fullfile(input.output_path, result_dir, "Input.mat"));
        all_output{id} = load(fullfile(input.output_path, result_dir, 'Output.mat'));
    end
    chain_summary(all_input, all_output, numData, input.output_path);
    fprintf('Completed chain search. \n')
end

function mustBeFileType(filename, ext)
    if ~endsWith(filename, ext)
        error("The provided name ('%s') must end with ('%s')", filename, ext)
    end
end
    
function data_name = getDataName(dataFolder)
    data_name = split(dataFolder, '\');
    
    if data_name(end) == ""
        data_name = data_name(1:end-1);
    end
    
    pattern = '^imec\d+_ks\d+$';
    if ~isempty(regexp(data_name(end), pattern, 'once'))
        data_name = data_name(end-1);
    else
        data_name = data_name(end);
    end
end
