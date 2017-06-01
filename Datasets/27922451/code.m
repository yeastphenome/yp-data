%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
vlaming_leeuwen_2016.pmid = 27922451;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(vlaming_leeuwen_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/elife-18919-fig2-data2-v2.xlsx', 'UpTag');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/elife-18919-fig2-data2-v2.xlsx', 'DownTag');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,3);
hit_strains2 = data2(2:end,3);

% Get the data itself
hit_data1 = cell2mat(data1(2:end,6));
hit_data2 = cell2mat(data2(2:end,6));

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);

hit_strains1(strncmp('BRE1PLATE', hit_strains1, length('BRE1PLATE'))) = {'YDL074C'};
hit_strains1(strncmp('DOT1PLATE', hit_strains1, length('DOT1PLATE'))) = {'YDR440W'};

hit_strains2(strncmp('BRE1PLATE', hit_strains2, length('BRE1PLATE'))) = {'YDL074C'};
hit_strains2(strncmp('DOT1PLATE', hit_strains2, length('DOT1PLATE'))) = {'YDR440W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

hit_strains1(inds) = [];
hit_data1(inds) = [];

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds)); 

hit_strains2(inds) = [];
hit_data2(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains),2);
[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1) = hit_data1(ind1);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,2) = hit_data2(ind1);

hit_data = nanmean(hit_data,2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5373];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
vlaming_leeuwen_2016.orfs = hit_strains;
vlaming_leeuwen_2016.ph = hit_data_names;
vlaming_leeuwen_2016.data = hit_data;
vlaming_leeuwen_2016.dataset_ids = hit_data_ids;

%% Save

save('./vlaming_leeuwen_2016.mat','vlaming_leeuwen_2016');

%% Print out

fid = fopen('./vlaming_leeuwen_2016.txt','w');
write_matrix_file(fid, vlaming_leeuwen_2016.orfs, vlaming_leeuwen_2016.ph, vlaming_leeuwen_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(vlaming_leeuwen_2016)
end

end

