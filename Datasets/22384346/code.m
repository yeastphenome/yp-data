%% Piggott~Nislow, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
piggott_nislow_2011.pmid = 22384346;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(piggott_nislow_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the hom data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS3.xlsx', 'ORF sort');

% Get the list of ORFs and the correponding data 
hom_strains = data(2:end,1);

% Get the data itself
hom_data = data(2:end,5:10); 
   
% Eliminate all white spaces & capitalize
hom_strains = clean_orf(hom_strains);

% If in gene name form, transform into ORF name
[hom_strains, translated, ambiguous] = translate(hom_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains));
hom_strains(inds) = [];
hom_data(inds, :) = [];

% Transform data
hom_data = cell2mat(hom_data);

% If the same strain is present more than once, average its values
[hom_strains, hom_data] = grpstats(hom_data, hom_strains, {'gname','mean'});

%% Load the het data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS6.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
het_strains = data(4:end,1);

% Get the data itself
het_data = data(4:end,4:9); 
   
% Eliminate all white spaces & capitalize
het_strains = clean_orf(het_strains);

% If in gene name form, transform into ORF name
[het_strains, translated, ambiguous] = translate(het_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(het_strains));
het_strains(inds) = [];
het_data(inds, :) = [];

% Transform data
het_data = cell2mat(het_data);

% If the same strain is present more than once, average its values
[het_strains, het_data] = grpstats(het_data, het_strains, {'gname','mean'});

%% Combine the datasets

% Add strains together to make final list of ORFs
hit_strains = unique([hom_strains; het_strains]);
hit_data = nan(length(hit_strains), 12);

% Intersect to make final data
[~, ind1, ind2] = intersect(hit_strains, hom_strains);
hit_data(ind1, 1:6) = hom_data(ind2, :);
[~, ind1, ind2] = intersect(hit_strains, het_strains);
hit_data(ind1, 7:12) = het_data(ind2, :);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11836; 11838; 11839; 11840; 11841; 11842; 11837; 11843; 11844; 11845; 11846; 11847];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
piggott_nislow_2011.orfs = hit_strains;
piggott_nislow_2011.ph = hit_data_names;
piggott_nislow_2011.data = hit_data;
piggott_nislow_2011.dataset_ids = hit_data_ids;

%% Save

save('./piggott_nislow_2011.mat','piggott_nislow_2011');

%% Print out

fid = fopen('./piggott_nislow_2011.txt','w');
write_matrix_file(fid, piggott_nislow_2011.orfs, piggott_nislow_2011.ph, piggott_nislow_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(piggott_nislow_2011)
end

end