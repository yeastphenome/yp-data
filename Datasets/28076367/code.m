%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bozaquel_morais_montero_lomeli_2017.pmid = 28076367;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bozaquel_morais_montero_lomeli_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pone.0169682.s002.xlsx', 'PRIMARY HITS');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(4:end,1);

% Get the data itself
hit_data = cell2mat(data(4:end,3:6));
hit_data = [sum(hit_data(:,1:2), 2), sum(hit_data(:,3:4),2)];

hit_data = -hit_data;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16260; 16261];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/Mat a collection YSC1053 open biosystems.xlsx', 'DATA');

tested_strains = tested_strains(2:end,2);

inds = find(cellfun(@isnumeric, tested_strains));
tested_strains(inds) = [];

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YLR287-A'})) = {'YLR287C-A'};

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
bozaquel_morais_montero_lomeli_2017.orfs = tested_strains;
bozaquel_morais_montero_lomeli_2017.ph = hit_data_names;
bozaquel_morais_montero_lomeli_2017.data = zeros(length(bozaquel_morais_montero_lomeli_2017.orfs),length(bozaquel_morais_montero_lomeli_2017.ph));
bozaquel_morais_montero_lomeli_2017.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, bozaquel_morais_montero_lomeli_2017.orfs);
bozaquel_morais_montero_lomeli_2017.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./bozaquel_morais_montero_lomeli_2017.mat','bozaquel_morais_montero_lomeli_2017');

%% Print out

fid = fopen('./bozaquel_morais_montero_lomeli_2017.txt','w');
write_matrix_file(fid, bozaquel_morais_montero_lomeli_2017.orfs, bozaquel_morais_montero_lomeli_2017.ph, bozaquel_morais_montero_lomeli_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bozaquel_morais_montero_lomeli_2017)
end

end

