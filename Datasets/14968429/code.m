%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
riles_reines_2004.pmid = 14968429;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(riles_reines_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/hits.txt','delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.ORF;

% Get the data itself
data.x6AUSensitivity(strcmp('Resistant', data.x6AUSensitivity)) = {'-1'};
hit_data = -cell2mat(cellfun(@str2num, data.x6AUSensitivity,'UniformOutput',0));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [188];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/aHOMDIP_FinalList_2804_15Aug03.xlsx', 'box order');

tested_strains = tested_strains(:,1);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YDL090'})) = {'YDL090C'};

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
riles_reines_2004.orfs = tested_strains;
riles_reines_2004.ph = hit_data_names;
riles_reines_2004.data = zeros(length(riles_reines_2004.orfs),length(riles_reines_2004.ph));
riles_reines_2004.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, riles_reines_2004.orfs);
riles_reines_2004.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./riles_reines_2004.mat','riles_reines_2004');

%% Print out

fid = fopen('./riles_reines_2004.txt','w');
write_matrix_file(fid, riles_reines_2004.orfs, riles_reines_2004.ph, riles_reines_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(riles_reines_2004)
end

end

