%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
krastel_hoepfner_2015.pmid = 26179970;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(krastel_hoepfner_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/Krastel_et_al_HIPHOPrawdata.xlsx', 'Nannocystin Exp 1');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/Krastel_et_al_HIPHOPrawdata.xlsx', 'Nannocystin Exp 2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,5);
hit_strains2 = data2(2:end,5);

type1 = data1(2:end,1);
type2 = data2(2:end,1);

% Get the data itself
% Keeping the data as is because, from the paper, it seems that negative
% values correspond to sensitivity (as per our convention)
hit_data1 = cell2mat(data1(2:end,3));
hit_data2 = cell2mat(data2(2:end,3));

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains),2,2);

inds = find(strcmp('HOP', type1));
[~,ind1,ind2] = intersect(hit_strains1(inds), hit_strains);
hit_data(ind2,1,1) = hit_data1(inds(ind1));
inds = find(strcmp('HIP', type1));
[~,ind1,ind2] = intersect(hit_strains1(inds), hit_strains);
hit_data(ind2,2,1) = hit_data1(inds(ind1));

inds = find(strcmp('HOP', type2));
[~,ind1,ind2] = intersect(hit_strains2(inds), hit_strains);
hit_data(ind2,1,2) = hit_data2(inds(ind1));
inds = find(strcmp('HIP', type2));
[~,ind1,ind2] = intersect(hit_strains2(inds), hit_strains);
hit_data(ind2,2,2) = hit_data2(inds(ind1));

% Average the replicates
hit_data = nanmean(hit_data,3);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5253 5254]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
krastel_hoepfner_2015.orfs = hit_strains;
krastel_hoepfner_2015.ph = hit_data_names;
krastel_hoepfner_2015.data = hit_data;
krastel_hoepfner_2015.dataset_ids = hit_data_ids;

%% Save

save('./krastel_hoepfner_2015.mat','krastel_hoepfner_2015');

%% Print out

fid = fopen('./krastel_hoepfner_2015.txt','w');
write_matrix_file(fid, krastel_hoepfner_2015.orfs, krastel_hoepfner_2015.ph, krastel_hoepfner_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(krastel_hoepfner_2015)
end

end

