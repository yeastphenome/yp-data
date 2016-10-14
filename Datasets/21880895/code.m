%% Shi~Emr, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
shi_emr_2011.pmid = 21880895;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(shi_emr_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_E11-05-0440_mc-E11-05-0440-s03.xlsx');

% Eliminate anything that doesn't look like an ORF
hit_orfs = data(4:end,2);

% Remove empty spaces and clean the orfs
hit_orfs(find(cellfun(@isnumeric, hit_orfs))) = [];
hit_orfs = clean_orf(hit_orfs);

% Remove arbitrary values
hit_orfs(strcmp('YPL072WA', hit_orfs)) = {'YPL072W-A'};
hit_orfs = unique(hit_orfs);

inds = find(~is_orf(hit_orfs));
hit_orfs(inds) = [];

% Make a data matrix of 1s
hit_data = zeros(size(hit_orfs))+1;

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [148];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

shi_emr_2011.orfs = hit_orfs;
shi_emr_2011.data = hit_data;
shi_emr_2011.ph = hit_data_names;
shi_emr_2011.dataset_ids = hit_data_ids;

%% Save

save('./shi_emr_2011.mat','shi_emr_2011');

%% Print out

fid = fopen('./shi_emr_2011.txt','w');
write_matrix_file(fid, shi_emr_2011.orfs, shi_emr_2011.ph, shi_emr_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(shi_emr_2011)
end

end
