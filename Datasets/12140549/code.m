%% Giaever~Johnston, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
giaever_johnston_2002.pmid = 12140549;

phenotypes = {'phenotype1'; 'phenotype2'};
treatments = {'condition1'; 'condition2'};


%% Hit strains

% Load the phenotype mapping file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/phenotype_mapping.xlsx','Sheet1');
pm.datasetid0 = cell2mat(data(2:end,1));
pm.who = data(2:end,7);
pm.what = data(2:end,8);
pm.when = data(2:end,6);
pm.how = data(2:end,9);
pm.condition = data(2:end,3);

for i = 1 : length(pm.condition)
    inds = regexp(pm.condition{i}, ' ');
    pm.condition{i} = pm.condition{i}(1:inds(end-2)-1);
end

pm = build_phenotype_name(pm);

[phenotypes, ia, ic] = unique(pm.ph);

hit_data_files = dir('./raw_data/');
hit_data_files_names = {hit_data_files(:).name}';
inds = find(cellfun(@isempty, regexp(hit_data_files_names,'^[0-9]{1,2}\_')));
hit_data_files_names(inds) = [];

phenotypes_hit_strains = cell(size(phenotypes));
phenotypes_hit_data = cell(size(phenotypes));

for i = 1 : length(phenotypes)
    
        % Find dataset_ids
        inds = find(ic == i);
        
        % Find the files
        e = strjoin(pm.datasetid0(inds),'|');
        indfiles = find(~cellfun(@isempty, regexp(hit_data_files_names, ['^(' e ')\_'])));

        all_hit_strains = [];
        all_hit_data = [];
        
        for j = 1 : length(indfiles)
            % Load hit strains
            [FILENAMES{end+1}, data] = read_data('textscan',['./raw_data/' hit_data_files_names{indfiles(j)}], '%s %f %*[^\n]');

            % Get the list of ORFs and the correponding data 
            hit_strains = data{1};

            % Get the data itself
            hit_data = exp(-data{2});   % this step is necessary to combine sensitivity and resistance data together (otherwise, they both have positive & negative values)
            if ~isempty(regexp(hit_data_files_names{indfiles(j)},'sen'))
                hit_data = -hit_data;
            end
    
            % Eliminate all white spaces & capitalize
            hit_strains = clean_orf(hit_strains);

            % Find anything that doesn't look like an ORF
            inds = find(~is_orf(hit_strains));
            disp(hit_strains(inds)); 
            
            all_hit_strains = [all_hit_strains; hit_strains];
            all_hit_data = [all_hit_data; hit_data];
        end

    % If the same strain is present more than once, average its values
    [all_hit_strains, all_hit_data] = grpstats(all_hit_data, all_hit_strains, {'gname','mean'});
    
    phenotypes_hit_strains{i} = all_hit_strains;
    phenotypes_hit_data{i} = all_hit_data;
end

% Combine all data together into a single matrix
all_orfs = unique(vertcat(phenotypes_hit_strains{:}));
all_data = nan(length(all_orfs),length(phenotypes));
for i = 1 : length(phenotypes)
    [~,ind1,ind2] = intersect(all_orfs, phenotypes_hit_strains{i});
    all_data(ind1,i) = phenotypes_hit_data{i}(ind2);
end

%% Prepare final dataset

% If the dataset is quantitative:
giaever_johnston_2002.orfs = all_orfs;
giaever_johnston_2002.ph = phenotypes;
giaever_johnston_2002.data = all_data;

%% Save

save('./giaever_johnston_2002.mat','giaever_johnston_2002');

end

