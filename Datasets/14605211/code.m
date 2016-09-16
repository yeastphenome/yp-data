%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ran_lau_2003.pmid = 14605211;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ran_lau_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('textscan','./raw_data/sensitive_hits.txt', '%s');
[FILENAMES{end+1}, data2] = read_data('textscan','./raw_data/resistant_hits.txt', '%s');

% Get the list of ORFs and the correponding data 
hit_strains = [data1(:,1); data2(:,1)];

% Get the data itself
hit_data = [-ones(size(data1(:,1))); ones(size(data2(:,1)))];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);
hit_strains(find(strcmp('FYV3', hit_strains))) = {'YDL151C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [140];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ran_lau_2003.orfs = hit_strains;
ran_lau_2003.ph = hit_data_names;
ran_lau_2003.data = hit_data;
ran_lau_2003.dataset_ids = hit_data_ids;

%% Save

save('./ran_lau_2003.mat','ran_lau_2003');

%% Print out

fid = fopen('./ran_lau_2003.txt','w');
write_matrix_file(fid, ran_lau_2003.orfs, ran_lau_2003.ph, ran_lau_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ran_lau_2003)
end

end

