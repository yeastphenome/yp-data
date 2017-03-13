%% Dunn~Jensen, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dunn_jensen_2006.pmid = 16267274;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(dunn_jensen_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Supp_Table3.xlsx', 'Sheet1');

ind_orf = strmatch('Systemic name', data.raw(12,:));
hit_orfs = upper(data.raw(14:end, ind_orf));

ind_data1 = strmatch('(-EtBr/+EtBr)', data.raw(12,:));
hit_data = 1./cell2mat(data.raw(14:end, ind_data1)); % reverse the ratio so that the lower the value, the sicker the mutant.

% Eliminate all white spaces & capitalize
hit_orfs = clean_orf(hit_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));  

% Average data for identical ORFs that appear multiple times
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [54];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
dunn_jensen_2006.orfs = hit_orfs;
dunn_jensen_2006.ph = hit_data_names;
dunn_jensen_2006.data = hit_data;
dunn_jensen_2006.dataset_ids = hit_data_ids;

%% Save

save('./dunn_jensen_2006.mat','dunn_jensen_2006');

%% Print out

fid = fopen('./dunn_jensen_2006.txt','w');
write_matrix_file(fid, dunn_jensen_2006.orfs, dunn_jensen_2006.ph, dunn_jensen_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(dunn_jensen_2006)
end

end
