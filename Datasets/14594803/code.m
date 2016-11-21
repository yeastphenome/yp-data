%% Davis-Kaplan~Kaplan, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
davis_kaplan_kaplan_2004.pmid = 14594803;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(davis_kaplan_kaplan_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

filename = './raw_data/Copy of Kaplan Experimental diploid homozygous deletion srceen Release 1 and 2 sort sheet(3).xlsx';
[~,sheets] = xlsfinfo(filename);
[FILENAMES{end+1}, data] = read_data('xlsread',filename, sheets{4});

% Get the list of ORFs and the correponding data 
hit_strains = data(3:end,3);

% Get the data itself
hit_data = data(3:end, 7:end); 
hit_screens = data(1,7:end)';

inds = find(cellfun(@isnumeric, hit_screens));
hit_screens(inds) = [];
hit_data(:,inds) = [];

hit_screens = strtrim(hit_screens);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains); 

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'TAL004W'})) = {'YAL004W'};
hit_strains(ismember(hit_strains, {'YELOO1C'})) = {'YEL001C'};
hit_strains(ismember(hit_strains, {'YCR102W-'})) = {'YCR102W'};
hit_strains(ismember(hit_strains, {'YDL133C-'})) = {'YDL133C'};
hit_strains(ismember(hit_strains, {'YKL096W-'})) = {'YKL096W'};
hit_strains(ismember(hit_strains, {'YAR002C-'})) = {'YAR002C'};
hit_strains(ismember(hit_strains, {'YBR084C-'})) = {'YBR084C'};
hit_strains(ismember(hit_strains, {'YCL026C-'})) = {'YCL026C'};
hit_strains(ismember(hit_strains, {'YCR028C-'})) = {'YCR028C'};
hit_strains(ismember(hit_strains, {'YOR298C-'})) = {'YOR298C'};
hit_strains(ismember(hit_strains, {'YPL183W-'})) = {'YPL183W'};
hit_strains(ismember(hit_strains, {'YHR132W-'})) = {'YHR132W'};
hit_strains(ismember(hit_strains, {'YKL053C-'})) = {'YKL053C'};
hit_strains(ismember(hit_strains, {'YOR298C-'})) = {'YOR298C'};
hit_strains(ismember(hit_strains, {'YFL013W-'})) = {'YFL013W'};
hit_strains(ismember(hit_strains, {'YLR390W-'})) = {'YLR390W'};
hit_strains(ismember(hit_strains, {'YFL013W-'})) = {'YFL013W'};
hit_strains(ismember(hit_strains, {'YLR390W-'})) = {'YLR390W'};
hit_strains(ismember(hit_strains, {'YBR090C-'})) = {'YBR090C'};
hit_strains(ismember(hit_strains, {'YBR162W-'})) = {'YBR162W'};
hit_strains(ismember(hit_strains, {'YER119C-'})) = {'YER119C'};
hit_strains(ismember(hit_strains, {'YFL010W-'})) = {'YFL010W'};
hit_strains(ismember(hit_strains, {'YPL183W-'})) = {'YPL183W'};
hit_strains(ismember(hit_strains, {'YGR122C-'})) = {'YGR122C'};
hit_strains(ismember(hit_strains, {'YPR133W-'})) = {'YPR133W'};
hit_strains(ismember(hit_strains, {'YIL015C-'})) = {'YIL015C'};
hit_strains(ismember(hit_strains, {'YDL045W-'})) = {'YDL045W'};
hit_strains(ismember(hit_strains, {'YER014C-'})) = {'YER014C'};
hit_strains(ismember(hit_strains, {'KL187C'})) = {'YKL187C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% Fix the data
hit_data(~cellfun(@isnumeric, hit_data)) = {NaN};
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out

[FILENAMES{end+1}, dt] = read_data('textscan','./extras/dataset_ids.txt','%s %d','delimiter','\t');
[Lia, Locb] = ismember(hit_screens, strtrim(dt{1}));

hit_data_ids = dt{2}(Locb);

[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
davis_kaplan_kaplan_2004.orfs = hit_strains;
davis_kaplan_kaplan_2004.ph = hit_data_names;
davis_kaplan_kaplan_2004.data = hit_data;
davis_kaplan_kaplan_2004.dataset_ids = hit_data_ids;

%% Save

save('./davis_kaplan_kaplan_2004.mat','davis_kaplan_kaplan_2004');

%% Print out

fid = fopen('./davis_kaplan_kaplan_2004.txt','w');
write_matrix_file(fid, davis_kaplan_kaplan_2004.orfs, davis_kaplan_kaplan_2004.ph, davis_kaplan_kaplan_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(davis_kaplan_kaplan_2004)
end

end
