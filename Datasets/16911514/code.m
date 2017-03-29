%% Kawahata~Iefuji, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kawahata_iefuji_2006.pmid = 16911514;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kawahata_iefuji_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

phenotypes = {'growth [spot assay]'};
treatments = {'lactic acid [5.1% w/v], pH [2.7]'; 'lactic acid [3.1% w/v], pH [2.9]'; 'HCl [0.28% w/v], pH [2.4]';'HCl [0.24% w/v], pH [2.6]'; 'acetic acid [0.5% w/v], pH [4.2]'; 'acetic acid [0.4% w/v], pH [4.3]'};

%% Load resistant
[FILENAMES{end+1}, hits_resistant.raw] = read_data('xlsread','./raw_data/hits.xlsx', 'Resistant');

hits_resistant_orfs = [hits_resistant.raw(:,1); hits_resistant.raw(:,6)];
hits_resistant_orfs = clean_orf(hits_resistant_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_resistant_orfs));
disp(hits_resistant_orfs(inds));  

hits_resistant_orfs(inds) = [];

hits_resistant_scores = zeros(length(hits_resistant_orfs),3);

% Scores
tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,3)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,8)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,1) = 1;

tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,4)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,9)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,2) = 1;

tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,5)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,10)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,3) = 1;

%% Load sensitive
[FILENAMES{end+1}, hits_sensitive.raw] = read_data('xlsread','./raw_data/hits.xlsx', 'Sensitive');

hits_sensitive_orfs = [hits_sensitive.raw(:,1); hits_sensitive.raw(:,6)];
hits_sensitive_orfs = clean_orf(hits_sensitive_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_sensitive_orfs));
disp(hits_sensitive_orfs(inds));  

hits_sensitive_orfs(inds) = [];

hits_sensitive_scores = zeros(length(hits_sensitive_orfs),3);

% Scores
tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,3)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,8)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,1) = -1;

tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,4)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,9)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,2) = -1;

tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,5)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,10)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,3) = -1;

% Check overlap between resistant and sensitive
length(intersect(hits_resistant_orfs(hits_resistant_scores(:,3)>0), hits_sensitive_orfs(hits_sensitive_scores(:,3)<0)));

%% Combine data

hit_strains = unique([hits_resistant_orfs; hits_sensitive_orfs]);
hit_data = zeros(length(hit_strains), 6);
[~,ind1,ind2] = intersect(hits_resistant_orfs, hit_strains);
hit_data(ind2,[1 3 5]) = hits_resistant_scores(ind1,:);
[~,ind1,ind2] = intersect(hits_sensitive_orfs, hit_strains);
hit_data(ind2,[2 4 6]) = hits_sensitive_scores(ind1,:);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [422 177 418 419 420 421]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kawahata_iefuji_2006.orfs = hit_strains;
kawahata_iefuji_2006.ph = hit_data_names;
kawahata_iefuji_2006.data = hit_data;
kawahata_iefuji_2006.dataset_ids = hit_data_ids;

%% Save

save('./kawahata_iefuji_2006.mat','kawahata_iefuji_2006');

%% Print out

fid = fopen('./kawahata_iefuji_2006.txt','w');
write_matrix_file(fid, kawahata_iefuji_2006.orfs, kawahata_iefuji_2006.ph, kawahata_iefuji_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kawahata_iefuji_2006)
end

end
