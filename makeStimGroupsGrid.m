function [target_groups]=makeStimGroupsGrid(grid_coordinates_planesXYZP,n_targets,saveDir, num_blocks, trials_per_block)

%% make stim groups

num_gid_locations = size(grid_coordinates_planesXYZP,1);

n_targets = n_targets; % targets per group

n_groups = trials_per_block*num_blocks; % was 840; % assuming a 1 sec iti

% pull random neurons but ensure that each neuron is targeted the same
% numer of times

targets_total = n_targets*n_groups;
shots_eachNeuron = floor(targets_total/num_gid_locations);
shots_eachNeuron_raw = (targets_total/num_gid_locations);

pull_from = [repmat(1:num_gid_locations,[1 shots_eachNeuron]), randsample(1:num_gid_locations,targets_total-(shots_eachNeuron*num_gid_locations))];

pw_distances_um = squareform(pdist(grid_coordinates_planesXYZP(:,1:2)))*3.7;


answer = questdlg(['Total ' num2str(num_gid_locations) ' targets, each hit ' num2str(shots_eachNeuron_raw) ' times. Continue?'], ...
    'YES', ...
    'NO');


% Handle response
switch answer
    case 'Yes'
        disp('making masks');
    case 'No'
       return
end


target_groups = [];
for group = 1:n_groups
    target_groups(group).group=group;
    
    distance_threshold = 500; % in um, radius
    
    success = 0;
    failure_counter = 0;
    while ~success
        seed_neuron = randsample(pull_from,1);
        
        distances = pw_distances_um(seed_neuron,:);
        within_range = find(distances<distance_threshold);
        available = intersect(pull_from,within_range);
        available(find(available==seed_neuron))=[];
        
        try
            other_members = randsample(available,n_targets-1);
            if numel(available) == 1;
                success = 0;
            else
                success = 1;
            end
        catch
            success = 0;
            % repeat pulling the seed neuron
            failure_counter = failure_counter+1;
        end
        
        if failure_counter > 100 % if no more gorups can be made based on the available neurons, relax distance threshold for the remaining few groups
            distance_threshold = distance_threshold+100;
        end
        if failure_counter > 500; % if after relaxing distance criteriion still no groups can be made, take whatever neurons are left (might make group size smaller)
            other_members = available;
            success = 1;
        end
    end
    
    all_members = [seed_neuron other_members];
    
    these_coordinatesXYZ = grid_coordinates_planesXYZP(all_members,1:3);
    [~, GroupCentroid] = ZOBlockAvoider(these_coordinatesXYZ(:,1:2));
    this_disp = mean(pdist2(these_coordinatesXYZ(:,1:2),GroupCentroid));
    
    target_groups(group).these=all_members';
    target_groups(group).these_coordinatesXYZ=these_coordinatesXYZ;
    target_groups(group).this_disp=this_disp;
    
    % delete pulled neurons from pool
    [~,ia]=intersect(pull_from,all_members);
    pull_from(ia) = [];
    
    disp(['group ' num2str(group) ' of ' num2str(n_groups) ' done'])
end
disp(['each grid location targeted ' num2str(shots_eachNeuron_raw) ' times'])
disp(['total: ' num2str(num_gid_locations) ' grid locations'])

save([saveDir filesep 'makeStimGroupsGrid_workspace' '.mat'])

end