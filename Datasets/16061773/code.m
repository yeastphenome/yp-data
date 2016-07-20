%% Hellauer~Turcotte, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hellauer_turcotte_2005.pmid = 16061773;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hellauer_turcotte_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

[FILENAMES{end+1}, DATA] = read_data('textread','./raw_data/hits_all.txt', '%s %d');

hits_orfs = DATA{1};
scores = DATA{2};

% Eliminate white spaces before/after ORF
hits_orfs = clean_orf(hits_orfs);
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

% Load tested genes
[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/ORF.txt', '%s');

% Eliminate anything that doesn't look like an ORF
tested_orfs = clean_orf(tested_orfs);
tested_orfs(strcmp('YPL072WA', tested_orfs)) = {'YPL072W-A'};
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));
tested_orfs(inds) = [];

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [184];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

hellauer_turcotte_2005.orfs = tested_orfs;
hellauer_turcotte_2005.data = zeros(length(tested_orfs), 1);
[~,ind1,ind2] = intersect(hellauer_turcotte_2005.orfs, hits_orfs);
hellauer_turcotte_2005.data(ind1,:) = scores(ind2,:);
hellauer_turcotte_2005.ph = hit_data_names;
hellauer_turcotte_2005.dataset_ids = hit_data_ids;

%% Save

save('./hellauer_turcotte_2005.mat','hellauer_turcotte_2005');

%% Print out

fid = fopen('./hellauer_turcotte_2005.txt','w');
write_matrix_file(fid, hellauer_turcotte_2005.orfs, hellauer_turcotte_2005.ph, hellauer_turcotte_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hellauer_turcotte_2005)
end

end
