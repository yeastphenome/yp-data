%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
fleming_blackman_2002.pmid = 11830665;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(fleming_blackman_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the sensitivity data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/hit_data.txt','delimiter','\t');
hit_strains = data.genenames;
hit_data = [data.ps_519 data.ps_341];

wt = repmat(hit_data(1,:), length(hit_strains),1);
hit_data = (hit_data-wt) ./ wt;
hit_strains(1) = [];
hit_data(1,:) = [];

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated] = translate(hit_strains);
hit_strains(~translated) = [];
hit_data(~translated,:) = [];

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Load the resistance data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/hit_data2.txt','delimiter','\t');

% Eliminate all white spaces & capitalize
hit_strains2 = clean_genename(data.genenames);

% If in gene name form, transform into ORF name
[hit_strains2, translated] = translate(hit_strains2);
hit_strains2(~translated) = [];

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

%%

hit_strains = [hit_strains; hit_strains2];
hit_data = [hit_data; ones(length(hit_strains2),2)];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1308 471]';

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/Table 2 (web supplement).xlsx', 'Sheet1');

% Eliminate all white spaces & capitalize
tested_strains = tested_strains(strcmp('+', tested_strains(:,3)),1);
tested_strains = clean_orf(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

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

% If the dataset is quantitative:
fleming_blackman_2002.orfs = tested_strains;
fleming_blackman_2002.ph = hit_data_names;
fleming_blackman_2002.data = zeros(length(tested_strains), length(hit_data_names));

[~,ind1,ind2] = intersect(tested_strains, hit_strains);
fleming_blackman_2002.data(ind1,:) = hit_data(ind2,:);

fleming_blackman_2002.dataset_ids = hit_data_ids;
%% Save

save('./fleming_blackman_2002.mat','fleming_blackman_2002');

%% Print out

fid = fopen('./fleming_blackman_2002.txt','w');
write_matrix_file(fid, fleming_blackman_2002.orfs, fleming_blackman_2002.ph, fleming_blackman_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(fleming_blackman_2002)
end

end

