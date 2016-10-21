%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
griffith_devine_2003.pmid = 12871900;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(griffith_devine_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hits.txt', '%s %d','delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data{1};

% Get the data itself
hit_data = data{2};
   
% Eliminate all white spaces & capitalize
tmp = regexp(hit_strains, ' ', 'split');
hit_strains = cellfun(@(x) x(1), tmp);
hit_strains = clean_genename(hit_strains);

hit_strains(strcmp('FYV3', hit_strains)) = {'BUD30'};
hit_strains(strcmp('KRE24', hit_strains)) = {'YPL102C'};
hit_strains(strcmp('SDF1', hit_strains)) = {'TIP41'};
hit_strains(strcmp('TCI1', hit_strains)) = {'ACL4'};


% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [480];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains1] = read_data('xlsread','./raw_data/Res Gen diploid knock01.xlsx', 'Res Gen diploid knock01.txt');
[FILENAMES{end+1}, tested_strains2] = read_data('xlsread','./raw_data/Res Gen diploid knockouts02.xlsx', 'Res Gen diploid knockouts02.txt');

tested_strains1_orfs = tested_strains1(3:end,2);
tested_strains1_plates = cell2mat(tested_strains1(3:end,5));

tested_strains2_orfs = tested_strains2(3:end,2);
tested_strains2_plates = cell2mat(tested_strains2(3:end,5));

inds1 = find(tested_strains1_plates >= 301 & tested_strains1_plates <= 349);
inds2 = find(tested_strains2_plates >= 301 & tested_strains2_plates <= 349);

tested_strains = unique([tested_strains1_orfs(inds1); tested_strains2_orfs(inds2)]);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'KL187C'})) = {'YKL187C'};
tested_strains(ismember(tested_strains, {'TAL004W'})) = {'YAL004W'};
tested_strains(ismember(tested_strains, {'YELOO1C'})) = {'YEL001C'};

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
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
griffith_devine_2003.orfs = tested_strains;
griffith_devine_2003.ph = hit_data_names;
griffith_devine_2003.data = zeros(length(griffith_devine_2003.orfs),length(griffith_devine_2003.ph));
griffith_devine_2003.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, griffith_devine_2003.orfs);
griffith_devine_2003.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./griffith_devine_2003.mat','griffith_devine_2003');

%% Print out

fid = fopen('./griffith_devine_2003.txt','w');
write_matrix_file(fid, griffith_devine_2003.orfs, griffith_devine_2003.ph, griffith_devine_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(griffith_devine_2003)
end

end

