%% Jaime~Nislow, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jaime_nislow_2012.pmid = 22727066;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jaime_nislow_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/12864_2012_4394_MOESM2_ESM.xlsx', 'Table_S1');

% Get the list of ORFs and the correponding data 
hit_strains = data(4:end,1);

% Split the names by semicolon
hit_strains = regexp(hit_strains, ':', 'split');
hit_strains = cellfun(@(x) x{1}, hit_strains, 'UniformOutput', false);

% Get the data itself
hit_data = data(4:end,2);
hit_data = -cell2mat(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Split hom and het data using the list of essential genes (no other way to
% do it for this dataset).
load essential_genes_151215.mat
hit_data2 = nan(length(hit_strains),2);

inds = find(~ismember(hit_strains, essential_genes));
hit_data2(inds,1) = hit_data(inds);

inds = find(ismember(hit_strains, essential_genes));
hit_data2(inds,2) = hit_data(inds);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5341 5342]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
jaime_nislow_2012.orfs = hit_strains;
jaime_nislow_2012.ph = hit_data_names;
jaime_nislow_2012.data = hit_data2;
jaime_nislow_2012.dataset_ids = hit_data_ids;

%% Save

save('./jaime_nislow_2012.mat','jaime_nislow_2012');

%% Print out

fid = fopen('./jaime_nislow_2012.txt','w');
write_matrix_file(fid, jaime_nislow_2012.orfs, jaime_nislow_2012.ph, jaime_nislow_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(jaime_nislow_2012)
end

end