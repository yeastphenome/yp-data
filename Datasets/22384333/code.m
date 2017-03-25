%% Hoon~Nislow, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hoon_nislow_2011.pmid = 22384333;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hoon_nislow_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_1.3.219_TableS3.xlsx', 'homozygous deletions');

% Get the list of ORFs
strains = data(5:end, 1);

% Clean up ORFs
strains = cellfun(@(x) strtok(x, ':'), strains, 'UniformOutput', false); 
strains = clean_orf(strains);

% Get data from hits
hit_data = -cell2mat(data(5:end, 6)); 

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
strains(inds) = [];
hit_data(inds) = [];

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [519; 755];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hoon_nislow_2011.orfs = strains;
hoon_nislow_2011.ph = hit_data_names;
hoon_nislow_2011.data = hit_data;
hoon_nislow_2011.dataset_ids = hit_data_ids;

%% Save

save('./hoon_nislow_2011.mat','hoon_nislow_2011');

%% Print out

fid = fopen('./hoon_nislow_2011.txt','w');
write_matrix_file(fid, hoon_nislow_2011.orfs, hoon_nislow_2011.ph, hoon_nislow_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hoon_nislow_2011)
end

end