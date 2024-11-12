% calculate L2, FR, and locations of all chains
function [L2_value, fr_all, fr_change, x_loc, z_loc] = chain_stats(all_input, all_output, metrics_name, full_chain, numData, output_path) %#ok<INUSD>
% all output needed for z_drift

% Preallocate arrays for efficiency
L2_value = NaN(size(full_chain,1), numData-1);
fr_all = NaN(size(full_chain,1), numData);
fr_fold_all = NaN(size(full_chain,1), numData);
fr_change = NaN(size(full_chain,1), numData-1);
z_loc = NaN(size(full_chain,1), numData);
x_loc = NaN(size(full_chain,1), numData);
% z_drift = NaN(size(full_chain,1), numData);

% loop through each chain
for ichain = 1:size(full_chain,1)
    for id = 1:numData-1
        
        % current clu
        clu_label1 = full_chain(ichain,id);
        clu_label2 = full_chain(ichain,id+1);
        
        pair_results = load(fullfile(all_input{id}.input.EMD_path,['EMD_post',num2str(id),'.mat']));
        % pair_results.all_results table has a line for each pair created
        % => number of lines = min( number in f1, number in f2).
        % column 2 = f2 label, column 3 = f1 label, column 4 = unweighted
        % L2 value.
        row_idx = find(pair_results.all_results(:,2) == clu_label2);
        L2_value(ichain,id) = pair_results.all_results(row_idx,4);
        
        
        % Calculate firing rate change
        % read excel summary
        tb1 = readtable(fullfile(all_input{id}.input.data_path1, metrics_name));
        tb2 = readtable(fullfile(all_input{id}.input.data_path2, metrics_name));
        
        % find FR
        fr1 = tb1.firing_rate;
        fr2 = tb2.firing_rate;
        fr_ave1 = sum(fr1)/size(tb1,1); %average across num clu
        fr_ave2 = sum(fr2)/size(tb2,1);
        
        
        % find all clusters
        clu1 = tb1.cluster_id+1;
        clu2 = tb2.cluster_id+1;
        
        % calculate fold change
        idx_row_ref = find(clu_label1 == clu1); % idx in all day 1 clusters
        fr_all(ichain,id) = fr1(idx_row_ref);
        fr_fold_all(ichain,id) = fr1(idx_row_ref)/fr_ave1; % how is unit is firing compare to average of this day
        idx_row_ref2 = find(clu_label2 == clu2);
        fr_all(ichain,id+1) = fr2(idx_row_ref2);
        fr_fold_all(ichain,id+1) = fr2(idx_row_ref2)/fr_ave2;
        fr_change(ichain,id) = (fr_fold_all(ichain,id+1)-fr_fold_all(ichain,id))/fr_fold_all(ichain,id); %percentage change of FR compare to prior day
        
        
        % find z locations
        % z_drift(ichain,id+1) = all_output{id}.output.z_mode;
        pair_results = load(fullfile(all_input{id}.input.EMD_path, all_input{id}.input.filename_post));
        f1 = pair_results.f1; % get all location data
        f2 = pair_results.f2;
        f1_labels = pair_results.f1_labels;
        f2_labels = pair_results.f2_labels;
        idx1 = find(f1_labels == clu_label1); % find clu index in all units
        idx2 = find(f2_labels == clu_label2);
        z_loc(ichain,id) = f1(idx1,2); % all z locations, col = days
        x_loc(ichain,id) = f1(idx1,1);
        if id == numData-1
            z_loc(ichain,id+1) = f2(idx2,2);
            x_loc(ichain,id+1) = f2(idx2,1);
        end
    end
end

save(fullfile(output_path,'chain_stats.mat'),'full_chain','L2_value','fr_all','fr_change','z_loc','x_loc')
end