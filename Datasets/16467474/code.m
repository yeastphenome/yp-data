%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chasse_dohlman_2006.pmid = 16467474;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(chasse_dohlman_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table_S1_copy.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = [data(2:end,2); data(2:end,5)];

% Get the data itself
hit_data = [data(2:end,3); data(2:end,6)];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Split the data into specific phenotypes
hit_data2 = zeros(length(hit_strains),3);
for i = 1 : length(hit_strains)
    tmp = regexp(hit_data{i},'[\s,]', 'split');
    for j = 2 : length(tmp)
        switch tmp{j}
            case 'I'
                hit_data2(i,1) = -0.5;
            case 'II'
                hit_data2(i,1) = 0.5;
            case 'III'
                hit_data2(i,1:2) = -1;
            case 'IV'
                hit_data2(i,2) = -0.5;
            case 'V'
                hit_data2(i,2) = 0.5;
            case 'VI'
                hit_data2(i,3) = 1;
        end
    end
end

hit_data = hit_data2;


% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11828 11829 5180]';

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/ Gene List_ResGenLibrary.xlsx', 'mat_a_101501.txt');

tested_strains = tested_strains(3:end, 3);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YLR287-A'})) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
chasse_dohlman_2006.orfs = tested_strains;
chasse_dohlman_2006.ph = hit_data_names;
chasse_dohlman_2006.data = zeros(length(chasse_dohlman_2006.orfs),length(chasse_dohlman_2006.ph));
chasse_dohlman_2006.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, chasse_dohlman_2006.orfs);
chasse_dohlman_2006.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./chasse_dohlman_2006.mat','chasse_dohlman_2006');

%% Print out

fid = fopen('./chasse_dohlman_2006.txt','w');
write_matrix_file(fid, chasse_dohlman_2006.orfs, chasse_dohlman_2006.ph, chasse_dohlman_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(chasse_dohlman_2006)
end

end

