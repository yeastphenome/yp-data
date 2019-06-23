%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bottoms_piotrowski_2018.pmid = 29329531;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bottoms_piotrowski_2018.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/YPDg+2.3%gval.v.YPD+gval.pw.rc.edger.2013_12_01_2 for spotfire 001.xlsx', 'YPDg+2.3%gval.v.YPD+gval.pw.rc.');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = cell2mat(data(2:end,3)); 
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Untranslated ORFs are all pseudogenes or blocked reading frames
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16211];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
bottoms_piotrowski_2018.orfs = hit_strains;
bottoms_piotrowski_2018.ph = hit_data_names;
bottoms_piotrowski_2018.data = hit_data;
bottoms_piotrowski_2018.dataset_ids = hit_data_ids;

%% Save

save('./bottoms_piotrowski_2018.mat','bottoms_piotrowski_2018');

%% Print out

fid = fopen('./bottoms_piotrowski_2018.txt','w');
write_matrix_file(fid, bottoms_piotrowski_2018.orfs, bottoms_piotrowski_2018.ph, bottoms_piotrowski_2018.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bottoms_piotrowski_2018)
end

end

