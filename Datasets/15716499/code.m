%% Deutschbauer~Giaever, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
deutschbauer_giaevere_2005.pmid = 15716499;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(deutschbauer_giaevere_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/OrfGeneData.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,[9 13 15 17]);
for i = 1:4
    temp = hit_data(:,i);
    temp(find(~cellfun(@isnumeric, temp))) = {NaN};
    hit_data(:,i) = temp;
end
hit_data = cell2mat(hit_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [521; 5257; 522; 5258];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
deutschbauer_giaevere_2005.orfs = hit_strains;
deutschbauer_giaevere_2005.ph = hit_data_names;
deutschbauer_giaevere_2005.data = hit_data;
deutschbauer_giaevere_2005.dataset_ids = hit_data_ids;

%% Save

save('./deutschbauer_giaevere_2005.mat','deutschbauer_giaevere_2005');

%% Print out

fid = fopen('./deutschbauer_giaevere_2005.txt','w');
write_matrix_file(fid, deutschbauer_giaevere_2005.orfs, deutschbauer_giaevere_2005.ph, deutschbauer_giaevere_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(deutschbauer_giaevere_2005)
end

end