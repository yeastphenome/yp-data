%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yazawa_uemura_2007.pmid = 17506111;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yazawa_uemura_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hits.txt', '%s');

hit_strains = data;

% Get the data itself
hit_data = -ones(size(hit_strains));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

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
hit_data_ids = [5261];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
yazawa_uemura_2007.orfs = hit_strains;
yazawa_uemura_2007.ph = hit_data_names;
yazawa_uemura_2007.data = hit_data;
yazawa_uemura_2007.dataset_ids = hit_data_ids;

%% Save

save('./yazawa_uemura_2007.mat','yazawa_uemura_2007');

%% Print out

fid = fopen('./yazawa_uemura_2007.txt','w');
write_matrix_file(fid, yazawa_uemura_2007.orfs, yazawa_uemura_2007.ph, yazawa_uemura_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yazawa_uemura_2007)
end

end

