%% Outten~Culotta, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
outten_culotta_2005.pmid = 15641941;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(outten_culotta_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/outten_culotta_2005.xlsx', 'Sheet1');

hits_orfs = data.raw(:,1);

% Eliminate white spaces before/after ORF
hits_orfs = clean_orf(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));  

hits_scores = cell2mat(data.raw(:,2));

% Adjust scores such that -1 = weakest phenotype, -4 = strongest phenotype.
hits_scores = hits_scores - 5;

% The 6 hits with no score (not sure why), set to 0
hits_scores(isnan(hits_scores)) = 0;

% If the same strain is present more than once, average its values
[hits_orfs, hits_scores] = grpstats(hits_scores, hits_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [110];

%% Load tested genes
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Yeast Knockout -BY4741.xlsx', 'mat_a_060701.txt');

tested_orfs = data.raw(2:end,2);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_orfs);

tested_strains(find(strcmp('YLR287-A', tested_strains))) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

tested_strains = unique(tested_strains);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_strains);

% From supplement: MED2 (YDL005C) not available in Mat-a background, was used Mat-alpha
% instead. So it should be added it to the list of tested
tested_strains = [tested_strains; {'YDL005C'}];

% YHR039C-A is not in the list of tested set, but it has the alias
% YHR039C-B, which is in the tested set. So, these are likely to be the same gene, so YHR039C-A should be renamed.
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
outten_culotta_2005.orfs = tested_strains;
outten_culotta_2005.ph = hit_data_names;
outten_culotta_2005.data = zeros(length(outten_culotta_2005.orfs),length(outten_culotta_2005.ph));
outten_culotta_2005.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, outten_culotta_2005.orfs);
outten_culotta_2005.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./outten_culotta_2005.mat','outten_culotta_2005');

%% Print out

fid = fopen('./outten_culotta_2005.txt','w');
write_matrix_file(fid, outten_culotta_2005.orfs, outten_culotta_2005.ph, outten_culotta_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(outten_culotta_2005)
end

end
