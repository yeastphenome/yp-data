%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jakubkova_tomaska_2016.pmid = 27711131;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jakubkova_tomaska_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pone.0164175.s006.xlsx', 'Tab S2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end,3);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds) = [];

hit_data2 = zeros(length(hit_strains),2);

inds = find(strcmp('Nig hypersensitive', hit_data));
hit_data2(inds,1) = -1;
inds = find(strcmp('Nig resistant', hit_data));
hit_data2(inds,1) = 1;
inds = find(strcmp('Val and Nig hypersensitive', hit_data));
hit_data2(inds,:) = -1;
inds = find(strcmp('Val and Nig resistant', hit_data));
hit_data2(inds,:) = 1;
inds = find(strcmp('Val hypersensitive', hit_data));
hit_data2(inds,2) = -1;
inds = find(strcmp('Val resistant', hit_data));
hit_data2(inds,2) = 1;

hit_data = hit_data2;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5182; 5176];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('readtable','./raw_data/strain_a_mating_type.txt', 'delimiter','\t');

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains.ORF_name);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
jakubkova_tomaska_2016.orfs = tested_strains;
jakubkova_tomaska_2016.ph = hit_data_names;
jakubkova_tomaska_2016.data = zeros(length(jakubkova_tomaska_2016.orfs),length(jakubkova_tomaska_2016.ph));
jakubkova_tomaska_2016.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, jakubkova_tomaska_2016.orfs);
jakubkova_tomaska_2016.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./jakubkova_tomaska_2016.mat','jakubkova_tomaska_2016');

%% Print out

fid = fopen('./jakubkova_tomaska_2016.txt','w');
write_matrix_file(fid, jakubkova_tomaska_2016.orfs, jakubkova_tomaska_2016.ph, jakubkova_tomaska_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(jakubkova_tomaska_2016)
end

end

