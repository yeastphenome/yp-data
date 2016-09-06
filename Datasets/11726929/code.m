%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bennett_resnick_2001.pmid = 11726929;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bennett_resnick_2001.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/sensitive_strains.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hit_strains = data(1:end,1);

% Get the data itself
hit_data = data(1:end,2);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data2 = nan(size(hit_data));

inds = find(strcmp('S/R', hit_data)); % as explain in the legend, "slow recovery but no loss in survival"
hit_data2(inds) = 0;

inds = find(~strcmp('S/R', hit_data));
hit_data2(inds) = -cellfun(@length, hit_data(inds));

hit_data = hit_data2;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = 465;

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
bennett_resnick_2001.orfs = hit_strains;
bennett_resnick_2001.ph = hit_data_names;
bennett_resnick_2001.data = hit_data;
bennett_resnick_2001.dataset_ids = hit_data_ids;

%% Save

save('./bennett_resnick_2001.mat','bennett_resnick_2001');

%% Print out

fid = fopen('./bennett_resnick_2001.txt','w');
write_matrix_file(fid, bennett_resnick_2001.orfs, bennett_resnick_2001.ph, bennett_resnick_2001.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bennett_resnick_2001)
end

end

