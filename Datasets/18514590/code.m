%% Serero~Boiteux, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
serero_boiteux_2008.pmid = 18514590;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(serero_boiteux_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

% Load tested

% ATTEMPT #1: go through the EXCEL files. PROBLEM: too many genes obtained
% (5617).
% home_dir = '/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Serero~Boiteux/Deletion_mutants_list/';
% folders = dir(home_dir);
% folders_names = {folders.name};
% folders_names(strncmp('.', folders_names,1)) = [];
%
% tested_orfs = [];
% for i = 1 : length(folders_names)
%     % Find all Excel files in the folder
%     xls_files = dir([home_dir folders_names{i} '/*.XLS']);
%     xls_files_names = {xls_files.name};
%
%     for j = 1 : length(xls_files_names)
%         filename = [home_dir folders_names{i} '/' xls_files_names{j}];
%         [status,sheets] = xlsfinfo(filename);
%
%         for k = 1 : length(sheets)
%             [tested.txt, tested.num, tested.raw] = xlsread(filename,sheets{k});
%             if ~isempty(tested.raw)
%                 inds = find(strcmp('BY4741', tested.raw(:,3)));
%                 tested_orfs = [tested_orfs; tested.raw(inds,2)];
%             end
%         end
%     end
%     i
% end

% ATTEMPT #2: GO through the DOC files.
% cd ~/Laboratory/Datasets/Phenotypes/2008_Serero~Boiteux/Deletion_mutants_list/
% Covert all DOC files into TXT files by running: sudo textutil -convert txt */*.DOC
% Read the TXT files that end with "1"


% Find all TXT files in the folder
txt_files = dir('./raw_data/*.txt');
txt_files_names = {txt_files.name};
rtf_files = dir('./raw_data/*.RTF');
txt_files_names = [txt_files_names {rtf_files.name}]';

tested_orfs = [];
for j = 1 : length(txt_files_names)
    % If filenames ends in "~1" or "a", load the list
    t = regexp(txt_files_names{j},'\.','split');
    if strcmp(t{1}(end),'1') | strcmp(t{1}(end),'a')
        [FILENAMES{end+1}, tst] = read_data('textread',['./raw_data/' txt_files_names{j}], '%s');
        tested_orfs = [tested_orfs; tst(is_orf(tst))];
    end
end

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs = unique(tested_orfs);

%% Load data
[FILENAMES{end+1}, DATA] = read_data('textread','./raw_data/hits_genenames.txt', '%s %s', 'delimiter', '\t');

hits_genenames = DATA{1};
hits_scores_txt = DATA{2};

hits_scores = cellfun(@length, hits_scores_txt)+1;
hits_scores = -hits_scores;

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);

[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];   % 25 ORFs added

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [99];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
serero_boiteux_2008.orfs = tested_orfs;
serero_boiteux_2008.ph = hit_data_names;
serero_boiteux_2008.data = zeros(length(serero_boiteux_2008.orfs),length(serero_boiteux_2008.ph));
serero_boiteux_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, serero_boiteux_2008.orfs);
serero_boiteux_2008.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./serero_boiteux_2008.mat','serero_boiteux_2008');

%% Print out

fid = fopen('./serero_boiteux_2008.txt','w');
write_matrix_file(fid, serero_boiteux_2008.orfs, serero_boiteux_2008.ph, serero_boiteux_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(serero_boiteux_2008)
end

end
