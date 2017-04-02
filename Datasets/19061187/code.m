%% Matsufuji~Nakagawa, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
matsufuji_nakagawa_2008.pmid = 19061187;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(matsufuji_nakagawa_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/S.c 5000.xlsx', 'remake');
tested_orfs = tested.raw(3:end,3);

tested_orfs = clean_orf(tested_orfs);

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = translate(tested_orfs);

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds))

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hits_orfs] = read_data('textread','./raw_data/hits_orfs.txt', '%s');

hits_orfs = clean_orfs(hits_orfs);

hits_orfs = translate(hits_orfs);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

% Correction: ORF YHR041C was in the hit list twice. One instance should
% have been YDR378C.

hits_orfs = [hits_orfs; {'YDR378C'}];
hits_orfs = unique(hits_orfs);

hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
tested_orfs = [tested_orfs; missing];   % 1 orfs to be added

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [165];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
matsufuji_nakagawa_2008.orfs = tested_orfs;
matsufuji_nakagawa_2008.ph = hit_data_names;
matsufuji_nakagawa_2008.data = zeros(length(matsufuji_nakagawa_2008.orfs),length(matsufuji_nakagawa_2008.ph));
matsufuji_nakagawa_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, matsufuji_nakagawa_2008.orfs);
matsufuji_nakagawa_2008.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./matsufuji_nakagawa_2008.mat','matsufuji_nakagawa_2008');

%% Print out

fid = fopen('./matsufuji_nakagawa_2008.txt','w');
write_matrix_file(fid, matsufuji_nakagawa_2008.orfs, matsufuji_nakagawa_2008.ph, matsufuji_nakagawa_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(matsufuji_nakagawa_2008)
end

end
