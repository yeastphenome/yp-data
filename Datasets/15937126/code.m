%% Sambade~Kane, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
sambade_kane_2005.pmid = 15937126;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(sambade_kane_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/ResGen 384 well set 14 plates.xlsx', 'ResGen MATa -384');
tested_orfs = tested.raw(4:end,2);

tested_orfs = clean_orf(tested_orfs);

tested_orfs(ismember(tested_orfs,{'YLR287-A'})) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds)); 

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

% Adjustments
hits_orfs(strcmpi('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];   % 2 orfs to be added

% If the same strain is present more than once, average its values
[hits_orfs, hits_scores] = grpstats(hits_scores, hits_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [179];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
sambade_kane_2005.orfs = tested_orfs;
sambade_kane_2005.ph = hit_data_names;
sambade_kane_2005.data = zeros(length(sambade_kane_2005.orfs),length(sambade_kane_2005.ph));
sambade_kane_2005.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, sambade_kane_2005.orfs);
sambade_kane_2005.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./sambade_kane_2005.mat','sambade_kane_2005');

%% Print out

fid = fopen('./sambade_kane_2005.txt','w');
write_matrix_file(fid, sambade_kane_2005.orfs, sambade_kane_2005.ph, sambade_kane_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(sambade_kane_2005)
end

end
