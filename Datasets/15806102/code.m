%% Giorgini~Muchowski, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
giorgini_muchowski_2005.pmid = 15806102;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(giorgini_muchowski_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, hits] = read_data('textread','./raw_data/giorgini_muchowski_2005_hits.txt', '%s');

% Eliminate white spaces before/after gene names
hits(:,1) = clean_genename(hits);

% Translate genenames to ORF
hits_orfs = translate(hits);

% Hits = LOF suppressors of Htt103Q sensitivity, so positive score.
scores = ones(length(hits_orfs),1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [185];

%% Load tested genes

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v2(1).0.xlsx', 'DATA');

tested_orfs = clean_orf(tested.raw(:,2));
tested_orfs(strcmp(tested_orfs,'YLR287-A')) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
giorgini_muchowski_2005.orfs = tested_orfs;
giorgini_muchowski_2005.ph = hit_data_names;
giorgini_muchowski_2005.data = zeros(length(giorgini_muchowski_2005.orfs),length(giorgini_muchowski_2005.ph));
giorgini_muchowski_2005.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, giorgini_muchowski_2005.orfs);
giorgini_muchowski_2005.data(ind2,:) = scores(ind1,:);

%% Save

save('./giorgini_muchowski_2005.mat','giorgini_muchowski_2005');

%% Print out

fid = fopen('./giorgini_muchowski_2005.txt','w');
write_matrix_file(fid, giorgini_muchowski_2005.orfs, giorgini_muchowski_2005.ph, giorgini_muchowski_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(giorgini_muchowski_2005)
end

end
