%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
deranieh_greenberg_2015.pmid = 26324718;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(deranieh_greenberg_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hit_list.txt', '%s %s %s','delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data{3};

hit_list2 = [];
for i = 1 : length(hit_strains)
    tmp = regexp(hit_strains{i}, ',', 'split');
    hit_list2 = [hit_list2; strtrim(tmp)'];
end

hit_list2(cellfun(@isempty, hit_list2)) = [];

hit_list3 = hit_list2;
for i = 1 : length(hit_list2)
    tmp = regexp(hit_list2{i}, '\/', 'split');
    hit_list3{i} = tmp{1};
end

hit_strains = hit_list3;

% Get the data itself
hit_data = -ones(size(hit_strains)); % if the dataset is binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

hit_strains(find(strcmp('HXT12', hit_strains))) = {'YIL170W'};
hit_strains(find(strcmp('SDL1', hit_strains))) = {'YIL167W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4955];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
deranieh_greenberg_2015.orfs = hit_strains;
deranieh_greenberg_2015.ph = hit_data_names;
deranieh_greenberg_2015.data = hit_data;
deranieh_greenberg_2015.dataset_ids = hit_data_ids;

%% Save

save('./deranieh_greenberg_2015.mat','deranieh_greenberg_2015');

%% Print out

fid = fopen('./deranieh_greenberg_2015.txt','w');
write_matrix_file(fid, deranieh_greenberg_2015.orfs, deranieh_greenberg_2015.ph, deranieh_greenberg_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(deranieh_greenberg_2015)
end

end

