%% Fei~Yang, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
fei_yang_2011.pmid = 21829381;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(fei_yang_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/hits.txt');

% Get the list of ORFs and the correponding data 
hit_strains = data.labels_row(:,1);

% Get the data itself
hit_data = data.data; 

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Normalize to WT
indWt = find(strcmp('WT', hit_strains));
hit_data = hit_data ./ repmat(hit_data(indWt,:),length(hit_strains),1);

hit_strains(indWt) = [];
hit_data(indWt,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [124; 700];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

fei_yang_2011.orfs = hit_strains;
fei_yang_2011.ph = hit_data_names;
fei_yang_2011.data = hit_data;
fei_yang_2011.dataset_ids = hit_data_ids;

%% Save

save('./fei_yang_2011.mat','fei_yang_2011');

fid = fopen('./fei_yang_2011.txt','w');
write_matrix_file(fid, fei_yang_2011.orfs, fei_yang_2011.ph, fei_yang_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(fei_yang_2011)
end

end

