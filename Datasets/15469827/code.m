%% Begley~Samson, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
begley_samson_2004.pmid = 15469827;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(begley_samson_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

phenotypes = {'growth (spot assay)'};
treatments = {'MMS';'t-BuOOH';'4NQO';'UV score'};

%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Begley2003.xlsx', 'Sheet1');

hits_orfs = data.raw(2:end,1);
hits_data = data.raw(2:end,7:10);

hits_orfs = clean_orf(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));  

% If possible, fix the problem (typos, omissions etc.)
hits_orfs(ismember(hits_orfs, {'YAR002AW'})) = {'YAR002W-A'};
hits_orfs(ismember(hits_orfs, {'YFL033AC'})) = {'YFL033C-A'};
hits_orfs(ismember(hits_orfs, {'YKLO72W'})) = {'YKL072W'};
hits_orfs(ismember(hits_orfs, {'YOLO57W'})) = {'YOL057W'};
hits_orfs(ismember(hits_orfs, {'YOLO62C'})) = {'YOL062C'};
hits_orfs(ismember(hits_orfs, {'YJL206-A'})) = {'YJL205C-A'};   
hits_orfs(ismember(hits_orfs, {'YLR287-A'})) = {'YLR287C-A'};   

hits_data = -cell2mat(hits_data);   % The sensitivity scores are such that the higher the number, the more sensitive the mutant.

[hits_orfs,hits_data] = grpstats(hits_data,hits_orfs,{'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [111 543 544 545]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
begley_samson_2004.orfs = hits_orfs;
begley_samson_2004.ph = hit_data_names;
begley_samson_2004.data = hits_data;
begley_samson_2004.dataset_ids = hit_data_ids;

%% Save

save('./begley_samson_2004.mat','begley_samson_2004');

%% Print out

fid = fopen('./begley_samson_2004.txt','w');
write_matrix_file(fid, begley_samson_2004.orfs, begley_samson_2004.ph, begley_samson_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(begley_samson_2004)
end

end
