%% Serero~Boiteux, 2008
function FILENAMES = code()
FILENAMES = {};

serero_boiteux_2008.pmid = 18514590;

phenotypes = {'growth (colony size)'};
treatments = {'CdCl2 [100 uM]'};

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

home_dir = './raw_data/Deletion_mutants_list/';
folders = dir(home_dir);
folders_names = {folders.name};
folders_names(strncmp('.', folders_names,1)) = [];

tested_orfs = [];
for i = 1 : length(folders_names)
% Find all TXT files in the folder
txt_files = dir([home_dir folders_names{i} '/*.txt']);
txt_files_names = {txt_files.name};
rtf_files = dir([home_dir folders_names{i} '/*.RTF']);
txt_files_names = [txt_files_names; {rtf_files.name}];

for j = 1 : length(txt_files_names)
% If filenames ends in "~1" or "a", load the list
t = regexp(txt_files_names{j},'\.','split');
if strcmp(t{1}(end),'1') | strcmp(t{1}(end),'a')
[FILENAMES{end+1}, tst] = dataread('textread',[home_dir folders_names{i} '/' txt_files_names{j}], '%s');
inds = find(strncmp('Y', tst,1));
tested_orfs = [tested_orfs; tst(inds)];
end
end
i;
end


inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
tested_orfs = upper(tested_orfs);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(tested_orfs);

% Load data
[FILENAMES{end+1}, DATA] = dataread('textread','./raw_data/hits_genenames.txt', '%s %s', 'delimiter', '\t');

hits_genenames = DATA{1};
hits_scores_txt = DATA{2};

hits_scores = cellfun(@length, hits_scores_txt)+1;
hits_scores = -hits_scores;

hits_genenames = upper(hits_genenames);
hits_orfs = genename2orf(hits_genenames,'noannot');

% Adjustments
hits_orfs(strcmpi('lys7', hits_orfs)) = {'YMR038C'};
hits_orfs(strcmpi('trf4', hits_orfs)) = {'YOL115W'};

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);


serero_boiteux_2008.orfs = tested_orfs;
serero_boiteux_2008.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
serero_boiteux_2008.data(ind2,1) = hits_scores(ind1);

serero_boiteux_2008.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'serero_boiteux_2008.mat'],'serero_boiteux_2008');
return;

% Save data into database
dt = serero_boiteux_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

