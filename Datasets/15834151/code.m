%% Deng~Brown, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
deng_brown_2005.pmid = 15834151;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(deng_brown_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES{end+1}, data] = read_data('readtable', './raw_data/table3.tab.txt', 'Delimiter', '\t');

hit_strains = data.ORF;
hit_data = data.T_C;

% Clean Orfs
hit_strains = clean_orf(hit_strains);

inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

% Average any repeats
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Prepare Final Dataset

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [581];

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
deng_brown_2005.orfs = hit_strains;
deng_brown_2005.ph = hit_data_names;
deng_brown_2005.data = hit_data;
deng_brown_2005.dataset_ids = hit_data_ids;

%% Save

save('./deng_brown_2005.mat','deng_brown_2005');

%% Print out

fid = fopen('./deng_brown_2005.txt','w');
write_matrix_file(fid, deng_brown_2005.orfs, deng_brown_2005.ph, deng_brown_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(deng_brown_2005)
end

end
