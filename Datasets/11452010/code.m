%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ni_snyder_2001.pmid = 11452010;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ni_snyder_2001.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textread','./raw_data/phenotypes.txt', '%s %s','delimiter','\t');

% Get the list of ORFs and the correponding data 
hit_strains = data{1};

% Get the data itself
hit_data = data{2};
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains2 = unique(hit_strains);
hit_phenotypes2 = unique(hit_data);
hit_data2 = zeros(length(hit_strains2),3);

for i = 1 : length(hit_phenotypes2)
    inds = find(strcmp(hit_phenotypes2{i}, hit_data));
    [~,ind1,ind2] = intersect(hit_strains(inds), hit_strains2);
    hit_data2(ind2,i) = 1;
end


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4921 4922 467]';

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains1] = read_data('xlsread','./raw_data/homo1.xls', 'Sheet1');
[FILENAMES{end+1}, tested_strains2] = read_data('xlsread','./raw_data/homo2.xls', 'Sheet1');
[FILENAMES{end+1}, tested_strains3] = read_data('xlsread','./raw_data/homo3.xls', 'Sheet1');
[FILENAMES{end+1}, tested_strains4] = read_data('xlsread','./raw_data/homo4.xls', 'Sheet1');

tested_strains = [tested_strains1(:,2); tested_strains2(:,2); tested_strains3(:,2); tested_strains4(:,2)];

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If possible, fix the typo
tested_strains(strcmp(tested_strains, 'TAL004W')) = {'YAL004W'};
tested_strains(strcmp(tested_strains, 'YELOO1C')) = {'YEL001C'};
tested_strains(strcmp(tested_strains, 'KL187C')) = {'YKL187C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

% If it seems reasonable, add the missing hits to the list of tested strains
tested_strains = [tested_strains; missing]; % 2 added strains

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
ni_snyder_2001.orfs = tested_strains;
ni_snyder_2001.ph = hit_data_names;
ni_snyder_2001.data = zeros(length(ni_snyder_2001.orfs),length(ni_snyder_2001.ph));
ni_snyder_2001.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains2, ni_snyder_2001.orfs);
ni_snyder_2001.data(ind2,:) = hit_data2(ind1,:);

%% Save

save('./ni_snyder_2001.mat','ni_snyder_2001');

%% Print out

fid = fopen('./ni_snyder_2001.txt','w');
write_matrix_file(fid, ni_snyder_2001.orfs, ni_snyder_2001.ph, ni_snyder_2001.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ni_snyder_2001)
end

end

