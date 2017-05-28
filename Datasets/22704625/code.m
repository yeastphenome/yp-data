%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hoepfner_winzeler_2012.pmid = 22704625;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hoepfner_winzeler_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/Hoepfner_et_al_HIPHOPrawdate.xlsx', 'Cladosporin Exp 1');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/Hoepfner_et_al_HIPHOPrawdate.xlsx', 'Cladosporin Exp 2');
[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/Hoepfner_et_al_HIPHOPrawdate.xlsx', 'Cladosporin Exp 3');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,5);
hit_strains2 = data2(2:end,5);
hit_strains3 = data3(2:end,5);


% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);
hit_strains3 = clean_orf(hit_strains3);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);
hit_strains3 = translate(hit_strains3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

inds = find(~is_orf(hit_strains3));
disp(hit_strains3(inds));  

% Get the data itself
hit_data1 = -cell2mat(data1(2:end,3));
hit_data2 = -cell2mat(data2(2:end,3));
hit_data3 = -cell2mat(data3(2:end,3));

hit_strains = unique([hit_strains1; hit_strains2; hit_strains3]);
hit_data = nan(length(hit_strains), 2, 3);

inds_hop = find(strcmp('HOP', data1(2:end,1)));
inds_hip = find(strcmp('HIP', data1(2:end,1)));

[~,ind1,ind2] = intersect(hit_strains1(inds_hop), hit_strains);
hit_data(ind2,1,1) = hit_data1(inds_hop(ind1),1);
[~,ind1,ind2] = intersect(hit_strains1(inds_hip), hit_strains);
hit_data(ind2,2,1) = hit_data1(inds_hip(ind1),1);

inds_hop = find(strcmp('HOP', data2(2:end,1)));
inds_hip = find(strcmp('HIP', data2(2:end,1)));

[~,ind1,ind2] = intersect(hit_strains2(inds_hop), hit_strains);
hit_data(ind2,1,2) = hit_data2(inds_hop(ind1),1);
[~,ind1,ind2] = intersect(hit_strains2(inds_hip), hit_strains);
hit_data(ind2,2,2) = hit_data2(inds_hip(ind1),1);

inds_hop = find(strcmp('HOP', data3(2:end,1)));
inds_hip = find(strcmp('HIP', data3(2:end,1)));

[~,ind1,ind2] = intersect(hit_strains3(inds_hop), hit_strains);
hit_data(ind2,1,3) = hit_data3(inds_hop(ind1),1);
[~,ind1,ind2] = intersect(hit_strains3(inds_hip), hit_strains);
hit_data(ind2,2,3) = hit_data3(inds_hip(ind1),1);

% Average the replicates
hit_data = nanmean(hit_data,3);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5251; 5252];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hoepfner_winzeler_2012.orfs = hit_strains;
hoepfner_winzeler_2012.ph = hit_data_names;
hoepfner_winzeler_2012.data = hit_data;
hoepfner_winzeler_2012.dataset_ids = hit_data_ids;

%% Save

save('./hoepfner_winzeler_2012.mat','hoepfner_winzeler_2012');

%% Print out

fid = fopen('./hoepfner_winzeler_2012.txt','w');
write_matrix_file(fid, hoepfner_winzeler_2012.orfs, hoepfner_winzeler_2012.ph, hoepfner_winzeler_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hoepfner_winzeler_2012)
end

end

