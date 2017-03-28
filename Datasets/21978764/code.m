%% Svensson~Samson, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
svensson_samson_2011.pmid = 21978764;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(svensson_samson_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/1752-0509-5-157-s1.xlsx', '2. Gi50 and R2 all strains');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
strains = data.raw(:,2);

% Eliminate all white spaces & capitalize
strains = clean_orf(strains);

% Eliminate anything that doesn't look like an ORF
strains(ismember(strains, {'YAR002AW'})) = {'YAR002W'};
strains(ismember(strains, {'YOLO57W'})) = {'YOL057W'};
strains(ismember(strains, {'YKLO72W'})) = {'YKL072W'};
strains(ismember(strains, {'YJL206-A'})) = {'YJL206C'};
strains(ismember(strains, {'YLR287-A'})) = {'YLR287C-A'};
strains(ismember(strains, {'YFL033AC'})) = {'YFL033C'};
strains(ismember(strains, {'YOLO62C'})) = {'YOL062C'};

inds = find(~is_orf(strains));
strains(inds) = [];
data.raw(inds,:) = [];

% If in gene name form, transform into ORF name
strains = translate(strains);


% Separate deletions from DAMP strains. Personal communication from
% Peter Svensson: deletions are on plates 1-57, DAMPs are on plates
% 301-311

dels = 'ABCDEFGH';
plate = data.raw(:,1);
for i = 1 : length(dels)
    plate = strtok(plate, dels(i));
end
plate = cellfun(@str2num, plate);

inds = find(plate > 57);
strains(inds) = [];
data.raw(inds,:) = [];

% Isolate the data columns
raw_data = data.raw(:,4);

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, raw_data));
raw_data(inds) = {NaN};
raw_data = cell2mat(raw_data);

% Average data for identical ORFs that appear multiple times
[strains, raw_data] = grpstats(raw_data, strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [28];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
svensson_samson_2011.orfs = strains;
svensson_samson_2011.ph = hit_data_names;
svensson_samson_2011.data = raw_data;
svensson_samson_2011.dataset_ids = hit_data_ids;

%% Save

save('./svensson_samson_2011.mat','svensson_samson_2011');

%% Print out

fid = fopen('./svensson_samson_2011.txt','w');
write_matrix_file(fid, svensson_samson_2011.orfs, svensson_samson_2011.ph, svensson_samson_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(svensson_samson_2011)
end

end

