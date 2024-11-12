% This is the second function of the algorithm
% Input: Unit assignment from multiple datasets
% Output: Summary of chains within distance threshold
function [chain_all,z_loc,len] = chain_summary(all_input,all_output,numData,output_path)

% Summarize chains
numChain = 50; % starting size, expand as needed
z_loc = NaN(numChain,numData); % hold z locations
chain_all = NaN(numChain, numData); % hold chains

count = 0;

for ir = 1:length(all_output)-1
    EMD_path = all_input{ir}.input.EMD_path;
    pair1 = all_output{ir}.output.all_results_post; % KSgood pair label in day 1
    match1 = load(fullfile(EMD_path, ['EMD_post', num2str(ir), '.mat']));
    z1 = match1.f1(:,2); % z-location d1
    z2 = match1.f2(:,2); % z-location d2
    z_label1 = match1.f1_labels;
    z_label2 = match1.f2_labels;
    
    pair2 = all_output{ir+1}.output.all_results_post; % KSgood pair label in day 2
    match2 = load(fullfile(EMD_path, ['EMD_post', num2str(ir+1), '.mat']));
    z3 = match2.f2(:,2); % z-location in day 3
    z_label3 = match2.f2_labels;
    
    % Data Cleaning
    allPair1 = pair1(pair1(:,7) <= all_input{ir}.input.threshold, :); % only include pair above threshold
    allPair2 = pair2(pair2(:,7) <= all_input{ir}.input.threshold, :);
    
    % Find Chains
    [same_clu,~,~] = intersect(allPair1(:,2),allPair2(:,3),'stable'); % clusters preserved between two comparisons/3 datasets
    for ic = 1:length(same_clu) % append
        if any(chain_all(:,ir+1) == same_clu(ic))
            idx_align = find(chain_all(:,ir+1)==same_clu(ic)); % align with previous tracked clusters
            chain_all(idx_align,ir+2) = allPair2((allPair2(:,3) == same_clu(ic)),2); % cluster label in day x
            z_loc(idx_align,ir+2) = z3((z_label3 == chain_all(idx_align,ir+2)));
        else %create a half chain
            count = count + 1;
            % Expand chain_all and z_loc if count exceeds current size
            if count > size(chain_all, 1)
                chain_all = [chain_all; NaN(numChain, numData)];
                z_loc = [z_loc; NaN(numChain, numData)];
            end
            chain_all(count,ir+1) = same_clu(ic);
            chain_all(count,ir) = allPair1(allPair1(:,2) == same_clu(ic),3);
            chain_all(count,ir+2) = allPair2(allPair2(:,3) == same_clu(ic),2); %cluster label in day x
            z_loc(count,ir+1) = z2(z_label2 == chain_all(count,ir+1));
            z_loc(count,ir) = z1(z_label1 == chain_all(count,ir));
            z_loc(count,ir+2) = z3(z_label3 == chain_all(count,ir+2));
        end
    end
end

% Exclude rows that are entirely NaN (no chain)
chain_all(all(isnan(chain_all), 2), :) = [];
z_loc(all(isnan(z_loc), 2), :) = [];
len = sum(~isnan(chain_all), 2); % chain lengths

% Subtract 1 from cluster id's to match with original cluster id's (added 1 because MATLAB is 1-indexed)
chain_all = chain_all - 1;

if count > 0
    save(fullfile(output_path,'chain_summary.mat'), 'all_input', 'all_output', 'chain_all', 'z_loc', 'len')
else
    fprintf('No chains found.\n');
end

end