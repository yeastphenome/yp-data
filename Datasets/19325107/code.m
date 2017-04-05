%% Jonikas, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jonikas_schuldiner_2009.pmid = 19325107;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jonikas_schuldiner_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/NIHMS201195-supplement-st1.xlsx','Table_S1');

orfs = data.raw(3:end,1);
raw_data = data.raw(3:end,11:end);
headings = data.raw(1:2, 11:19);

orfs = clean_orf(orfs);

orfs(strcmp('YOLO62C', orfs)) = {'YOL062C'};
orfs(strcmp('YKLO72W', orfs)) = {'YKL072W'};
orfs(strcmp('YOLO57W', orfs)) = {'YOL057W'};
orfs(strcmp('YLR287-A', orfs)) = {'YLR287C-A'};

orfs = translate(orfs);

inds = find(~is_orf(orfs));
disp(orfs(inds));

orfs(inds) = [];
raw_data(inds,:) = [];

raw_data = raw_data(:,7);
headings = headings(:,7);

raw_data = cell2mat(raw_data);

[orfs, raw_data] = grpstats(raw_data, orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [699];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
jonikas_schuldiner_2009.orfs = orfs;
jonikas_schuldiner_2009.ph = hit_data_names;
jonikas_schuldiner_2009.data = raw_data;
jonikas_schuldiner_2009.dataset_ids = hit_data_ids;

%% Save

save('./jonikas_schuldiner_2009.mat','jonikas_schuldiner_2009');

%% Print out

fid = fopen('./jonikas_schuldiner_2009.txt','w');
write_matrix_file(fid, jonikas_schuldiner_2009.orfs, jonikas_schuldiner_2009.ph, jonikas_schuldiner_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(jonikas_schuldiner_2009)
end

end