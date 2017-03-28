%% Dos Santos~Sa-Correia, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

dos_santos_sa_correia_2011.pmid = 21960436;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(dos_santos_sa_correia_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the tested data

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx', 'Tabelle2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
tested_orfs = tested.raw(2:end,1);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% If in gene name form, transform into ORF name
tested_orfs = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];  

%% Load data
% Hypersensitive
[FILENAMES{end+1}, hits_genenames_HS] = read_data('textread','./raw_data/hits_genenames_hs.txt', '%s');

% Eliminate all white spaces & capitalize
hits_genenames_HS = clean_genename(hits_genenames_HS);

% If in gene name form, transform into ORF name
hits_orfs_HS = translate(hits_genenames_HS);

% Get the data itself
hits_scores_HS = zeros(length(hits_orfs_HS),1)-2;

% Sensitive
[FILENAMES{end+1}, hits_genenames_S] = read_data('textread','./raw_data/hits_genenames_s.txt', '%s');

% Eliminate all white spaces & capitalize
hits_genenames_S = clean_genename(hits_genenames_S);

% If in gene name form, transform into ORF name
hits_orfs_S = translate(hits_genenames_S);

% Get the data itself
hits_scores_S = zeros(length(hits_orfs_S),1)-1;

% Resistant
[FILENAMES{end+1}, hits_genenames_R] = read_data('textread','./raw_data/hits_genenames_r.txt', '%s');

% Eliminate all white spaces & capitalize
hits_genenames_R = upper(clean_genename(hits_genenames_R));

% If in gene name form, transform into ORF name
[hits_orfs_R, translated] = translate(hits_genenames_R);
hits_orfs_R(~translated) = [];

% Get the data itself
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

% Adjust the overlapping strains as follows: eliminate HS from S; eliminate
% HS ^ R and S ^ R from both
[~,~,ind2] = intersect(hits_orfs_HS, hits_orfs_S);
hits_orfs_S(ind2) = []; hits_scores_S(ind2) = [];

[~,ind1,ind2] = intersect(hits_orfs_HS, hits_orfs_R);
hits_orfs_HS(ind1) = []; hits_scores_HS(ind1) = [];
hits_orfs_R(ind2) = []; hits_scores_R(ind2) = [];

[~,ind1,ind2] = intersect(hits_orfs_S, hits_orfs_R);
hits_orfs_S(ind1) = []; hits_scores_S(ind1) = [];
hits_orfs_R(ind2) = []; hits_scores_R(ind2) = [];

hits_orfs = [hits_orfs_HS; hits_orfs_S; hits_orfs_R];
hits_scores = [hits_scores_HS; hits_scores_S; hits_scores_R];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Find anything that doesn't look like an ORF
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};
hits_orfs(strcmp('YBR089C-A', hits_orfs)) = {'YBR090C-A'};
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [149];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitive
dos_santos_sa_correia_2011.orfs = tested_orfs;
dos_santos_sa_correia_2011.ph = hit_data_names;
dos_santos_sa_correia_2011.data = zeros(length(tested_orfs), length(dos_santos_sa_correia_2011.ph));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
dos_santos_sa_correia_2011.data(ind2,1) = hits_scores(ind1);

%% Save

save('./dos_santos_sa_correia_2011.mat','dos_santos_sa_correia_2011');

%% Print out

fid = fopen('./dos_santos_sa_correia_2011.txt','w');
write_matrix_file(fid, dos_santos_sa_correia_2011.orfs, dos_santos_sa_correia_2011.ph, dos_santos_sa_correia_2011.data);
fclose(fid);


%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(firstauthor_lastauthor_YYYY)
end

end
