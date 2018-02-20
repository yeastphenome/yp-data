%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mattiazzi_petrovic_2010.pmid = 20379744;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mattiazzi_petrovic_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/ATX_SDL_combined_MattiazziUsaj_MGG2010.xlsx', 'ATX_SDL_combined_MattiazziUsaj_');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,5);

% Get the data itself
hit_data = data(2:end,8); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});


%% Load data 2 (discrete phenotypes)

[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/438_2010_533_MOESM2_ESM.xlsx', 'SDL screen 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2(2:end,1);

% Get the data itself
hit_data2 = data2(2:end,3); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

inds = find(cellfun(@isnumeric, hit_strains2));
hit_strains2(inds) = [];
hit_data2(inds) = [];

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_data2 = cell2mat(hit_data2);

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

%%

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16127;16126];

%% Merge the data

hit_strains_all = unique([hit_strains; hit_strains2]);
hit_data_all = nan(length(hit_strains_all), 2);

[~,ind1,ind2] = intersect(hit_strains, hit_strains_all);
hit_data_all(ind2,1) = hit_data(ind1);

[~,ind1,ind2] = intersect(hit_strains2, hit_strains_all);
hit_data_all(ind2,2) = hit_data2(ind1);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
mattiazzi_petrovic_2010.orfs = hit_strains_all;
mattiazzi_petrovic_2010.ph = hit_data_names;
mattiazzi_petrovic_2010.data = hit_data_all;
mattiazzi_petrovic_2010.dataset_ids = hit_data_ids;


%% Save

save('./mattiazzi_petrovic_2010.mat','mattiazzi_petrovic_2010');

%% Print out

fid = fopen('./mattiazzi_petrovic_2010.txt','w');
write_matrix_file(fid, mattiazzi_petrovic_2010.orfs, mattiazzi_petrovic_2010.ph, mattiazzi_petrovic_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mattiazzi_petrovic_2010)
end

end

