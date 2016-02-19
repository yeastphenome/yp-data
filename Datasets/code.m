%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
firstauthor_lastauthor_YYYY.pmid = 12345678;

phenotypes = {'phenotype1'; 'phenotype2'};
treatments = {'condition1'; 'condition2'};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/data.xlsx', 'Spreadsheet name');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = ones(size(hit_strains)); % if the dataset is binary
hit_data = data(:,2:4); % if the dataset is discrete or binary
    
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YAL001'})) = {'YAL001C'};

% If not possible, eliminate the entry
hit_strains(ismember(hit_strains, {'BLANK'})) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/tested_strains.xlsx', 'Spreadsheet name');

% Eliminate all white spaces & capitalize
tested_strains = clean_genename(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YAL001'})) = {'YAL001C'};

% If not possible, eliminate the entry
tested_strains(ismember(tested_strains, {'BLANK'})) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

% If it seems reasonable, add the missing hits to the list of tested strains
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% If the dataset is quantitative:
firstauthor_lastauthor_YYYY.orfs = hit_strains;
firstauthor_lastauthor_YYYY.ph = strcat(phenotypes, '; ', treatments);
firstauthor_lastauthor_YYYY.data = hit_data;


% If the dataset is discrete/binary and the tested strains were provided separately:
firstauthor_lastauthor_YYYY.orfs = tested_strains;
firstauthor_lastauthor_YYYY.ph = strcat(phenotypes, '; ', treatments);
firstauthor_lastauthor_YYYY.data = zeros(length(firstauthor_lastauthor_YYYY.orfs),length(firstauthor_lastauthor_YYYY.ph));

[~,ind1,ind2] = intersect(hit_strains, firstauthor_lastauthor_YYYY.orfs);
firstauthor_lastauthor.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./firstauthor_lastauthor_YYYY.mat','firstauthor_lastauthor_YYYY');

end

