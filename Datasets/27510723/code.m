%% Kwon~Koo, 2016
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kwon_koo_2016.pmid = 27510723;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kwon_koo_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the hom data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/425_2016_2579_MOESM3_ESM.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hom_strains = data(4:end,1);

% Get the data itself
hom_data = cell2mat(data(4:end,2));
   
% Transform the data
hom_data = hom_data .^ 10;
hom_data = 1 ./ hom_data;
hom_data = log(hom_data);

% Eliminate all white spaces & capitalize
hom_strains = clean_orf(hom_strains);

% If in gene name form, transform into ORF name
[hom_strains, translated, ambiguous] = translate(hom_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains));
hom_strains(inds) = [];
hom_data(inds) = [];

% If the same strain is present more than once, average its values
[hom_strains, hom_data] = grpstats(hom_data, hom_strains, {'gname','mean'});

%% Load the het data

[FILENAMES{end+1}, data] = read_data('xlsread',...
    './raw_data/het/Background Corrected Ratios And Robust Z-Scores (FD Scores)(all data).xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
het_strains = data(2:end,1);

% Get the data itself
het_data = data(2:end,17);
het_data(strcmp(het_data, {'NA'})) = {NaN};
het_data = cell2mat(het_data);
   
% Transform the data
het_data = het_data .^ 10;
het_data = 1 ./ het_data;
het_data = log(het_data);

% Eliminate all white spaces & capitalize
het_strains = clean_orf(het_strains);

% If in gene name form, transform into ORF name
[het_strains, translated, ambiguous] = translate(het_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(het_strains));
het_strains(inds) = [];
het_data(inds) = [];

% If the same strain is present more than once, average its values
[het_strains, het_data] = grpstats(het_data, het_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4826; 5181];

%% Combine the het and hom data
hit_strains = unique([het_strains; hom_strains]);
hit_data = nan(length(hit_strains),2);
[~, a, b] = intersect(hit_strains, het_strains);
hit_data(a, 2) = het_data(b);
[~, a, b] = intersect(hit_strains, hom_strains);
hit_data(a, 1) = hom_data(b);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kwon_koo_2016.orfs = hit_strains;
kwon_koo_2016.ph = hit_data_names;
kwon_koo_2016.data = hit_data;
kwon_koo_2016.dataset_ids = hit_data_ids;

%% Save

save('./kwon_koo_2016.mat','kwon_koo_2016');

%% Print out

fid = fopen('./kwon_koo_2016.txt','w');
write_matrix_file(fid, kwon_koo_2016.orfs, kwon_koo_2016.ph, kwon_koo_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kwon_koo_2016)
end

end