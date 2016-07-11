%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kubota_hirata_2004.pmid = 15118337;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kubota_hirata_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hits.txt', '%s');

hit_orfs = cell(size(data));
t = regexp(data, '/', 'split');
for i = 1 : length(t)
    hit_orfs{i} = t{i}{end};
end
   
% Eliminate all white spaces & capitalize
hit_orfs = clean_orf(hit_orfs);

% Find anything that doesn't look like an ORF
hit_orfs(strcmp('VHR060W', hit_orfs)) = {'YHR060W'};
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));  

hit_orfs = unique(hit_orfs);

hit_data = -ones(size(hit_orfs));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [191];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kubota_hirata_2004.orfs = hit_orfs;
kubota_hirata_2004.ph = hit_data_names;
kubota_hirata_2004.data = hit_data;
kubota_hirata_2004.dataset_ids = hit_data_ids;

%% Save

save('./kubota_hirata_2004.mat','kubota_hirata_2004');

%% Print out

fid = fopen('./kubota_hirata_2004.txt','w');
write_matrix_file(fid, kubota_hirata_2004.orfs, kubota_hirata_2004.ph, kubota_hirata_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kubota_hirata_2004)
end

end

