%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lis_bobek_2009.pmid = 19596888;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lis_bobek_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Lis_2009.xlsx', 'Log2 Ratio');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = cell2mat(data(2:end,2:4));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

hit_data = -hit_data; % log2 control/treatment must be reversed

% Normalize by 0 h
hit_data = hit_data - repmat(hit_data(:,1),1,3);
hit_data(:,1) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4924; 5255];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lis_bobek_2009.orfs = hit_strains;
lis_bobek_2009.ph = hit_data_names;
lis_bobek_2009.data = hit_data;
lis_bobek_2009.dataset_ids = hit_data_ids;

%% Save

save('./lis_bobek_2009.mat','lis_bobek_2009');

%% Print out

fid = fopen('./lis_bobek_2009.txt','w');
write_matrix_file(fid, lis_bobek_2009.orfs, lis_bobek_2009.ph, lis_bobek_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(lis_bobek_2009)
end

end

