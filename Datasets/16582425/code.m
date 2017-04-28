%% Hancock~Lopes, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hancock_lopes_2006.pmid = 16582425;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hancock_lopes_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/mat_alpha_061101.xlsx', 'mat_alpha_061101');

tested_orfs = tested.raw(4:end,2);
tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YYKL138C', tested_orfs)) = {'YKL138C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds)); 

tested_orfs(inds) = [];

tested_orfs = translate(tested_orfs);
tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hits_genenames] = read_data('textscan', './raw_data/hits_genenames.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
[hits_orfs, translated] = translate(hits_genenames);

hits_orfs = unique(hits_orfs);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hits_orfs, tested_orfs);
disp(missing);

hits_scores = ones(length(hits_orfs),1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [178];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
hancock_lopes_2006.orfs = tested_orfs;
hancock_lopes_2006.ph = hit_data_names;
hancock_lopes_2006.data = zeros(length(hancock_lopes_2006.orfs),length(hancock_lopes_2006.ph));
hancock_lopes_2006.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, hancock_lopes_2006.orfs);
hancock_lopes_2006.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./hancock_lopes_2006.mat','hancock_lopes_2006');

%% Print out

fid = fopen('./hancock_lopes_2006.txt','w');
write_matrix_file(fid, hancock_lopes_2006.orfs, hancock_lopes_2006.ph, hancock_lopes_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hancock_lopes_2006)
end

end
