%% Cooper~Fields, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
cooper_fields_2010.pmid = 20610602;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(cooper_fields_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/SupplementalTable4.xlsx');

% Get the list of ORFs and the correponding data 
hit_orfs = data.raw(2:end,1);

% Eliminate all white spaces & capitalize
hit_orfs = clean_orf(hit_orfs);

% Translate
[hit_orfs, translated, ambiguous]  = translate(hit_orfs);

% Find anything that doesn't look like an ORF
hit_orfs(strcmp('YML048WA-', hit_orfs)) = {'YML048W-A'};
inds = find(~is_orf(hit_orfs));
hit_orfs(inds) = [];

% Replace all void data values
inds = find(strcmp('-', data.raw));
data.raw(inds) = {NaN};

% Transform to numeric array
hit_data = cell2mat(data.raw(2:end,3:end));

% Combine lysine columns
lys = mean(hit_data(:,13:2:15), 2, 'omitnan'); 
hit_data(:,15) = [];
hit_data(:,13) = lys;

% Identical ORFs with multiple values -> average
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [32; 730; 731; 732; 733; 734; 735; 736; 737; 738; 739; 740; 741; 742; 743; 744];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
cooper_fields_2010.orfs = hit_orfs;
cooper_fields_2010.ph = hit_data_names;
cooper_fields_2010.data = hit_data;
cooper_fields_2010.dataset_ids = hit_data_ids;

%% Save

save('./cooper_fields_2010.mat','cooper_fields_2010');

%% Print out

fid = fopen('./cooper_fields_2010.txt','w');
write_matrix_file(fid, cooper_fields_2010.orfs, cooper_fields_2010.ph, cooper_fields_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(cooper_fields_2010)
end

end

