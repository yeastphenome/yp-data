%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
pan_boeke_2004.pmid = 15525520;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(pan_boeke_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS2.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(4:end,1);

% Get the data itself
hit_data = data(4:end,2:10); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = -cell2mat(hit_data);  % reversing the sign to indicate that lower numbers = slower growth

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5237:5245]';

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
% Microarray platform
geodata = getgeodata('GPL1444', 'ToFile','./raw_data/GPL1444.txt');
geodata.id = cell2mat(geodata.Data(:,1));
geodata.orf = geodata.Data(:,6);

% Untreated control
geodata2 = getgeodata('GSM30549', 'ToFile', './raw_data/GSM30549.txt');
geodata2.id = geodata2.Data(:,1);
geodata2.flag = geodata2.Data(:,9);

% Get all ORFs that have a flag == 0 (passed the basic filter)
geodata2.orf = cell(size(geodata2.id));
[~,ind1,ind2] = intersect(geodata.id, geodata2.id);
geodata2.orf(ind2) = geodata.orf(ind1);

inds = find(cellfun(@isnumeric, geodata2.orf));
geodata2.orf(inds) = [];
geodata2.id(inds) = [];
geodata2.flag(inds) = [];

geodata2.orf = clean_orf(geodata2.orf);
geodata2.orf2 = translate(geodata2.orf);

[orfs, flags] = grpstats(geodata2.flag, geodata2.orf2, {'gname','mean'});
inds = find(flags > -50);
tested_strains = orfs(inds);

% Excluded strains
[FILENAMES{end+1}, excluded] = read_data('xlsread','./raw_data/TableS3.xlsx', 'Table 1');
excluded_strains = excluded(2:end,1);

excluded_strains = clean_orf(excluded_strains);
excluded_strains = translate(excluded_strains);

inds = find(~is_orf(excluded_strains));
disp(excluded_strains(inds));

excluded_strains(inds) = [];
excluded_strains = unique(excluded_strains);

tested_strains = setdiff(tested_strains, excluded_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
pan_boeke_2004.orfs = tested_strains;
pan_boeke_2004.ph = hit_data_names;
pan_boeke_2004.data = zeros(length(pan_boeke_2004.orfs),length(pan_boeke_2004.ph));
pan_boeke_2004.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, pan_boeke_2004.orfs);
pan_boeke_2004.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./pan_boeke_2004.mat','pan_boeke_2004');

%% Print out

fid = fopen('./pan_boeke_2004.txt','w');
write_matrix_file(fid, pan_boeke_2004.orfs, pan_boeke_2004.ph, pan_boeke_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(pan_boeke_2004)
end

end

