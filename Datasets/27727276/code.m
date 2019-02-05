%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
payen_dunham_2016.pmid = 27727276;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(payen_dunham_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pgen.1006339.s008.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(10:end,1);

% Get the data itself
hit_data = cell2mat(data(10:end,2:14));
hit_data_names = data(9,2:14)';
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

datasets_to_keep = {'MM1N-Phosphate','MM2N-Phosphate','MM1N-Sulfate','MM2N-Sulfate','MM1N-Glucose','MM2N-Glucose'};
[~,ind1, ind2] = intersect(hit_data_names, datasets_to_keep);
[~,ix] = sort(ind2);
ind1 = ind1(ix);

hit_data = hit_data(:,ind1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16195; 16199; 16210; 16197; 16196; 16198];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
payen_dunham_2016.orfs = hit_strains;
payen_dunham_2016.ph = hit_data_names;
payen_dunham_2016.data = hit_data;
payen_dunham_2016.dataset_ids = hit_data_ids;

%% Save

save('./payen_dunham_2016.mat','payen_dunham_2016');

%% Print out

fid = fopen('./payen_dunham_2016.txt','w');
write_matrix_file(fid, payen_dunham_2016.orfs, payen_dunham_2016.ph, payen_dunham_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(payen_dunham_2016)
end

end

