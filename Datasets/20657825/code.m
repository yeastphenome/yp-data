%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
fabrizio_longo_2010.pmid = 20657825;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(fabrizio_longo_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pgen.1001024.s006.xlsx', 'SuppleTable2_v2.txt');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(4:end,1);

% Get the data itself
hit_data = cell2mat(data(4:end,5:12));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4711 4781 4782 4783 4711 4781 4782 4783]';
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
fabrizio_longo_2010.orfs = hit_strains;
fabrizio_longo_2010.ph = hit_data_names;
fabrizio_longo_2010.data = hit_data;
fabrizio_longo_2010.dataset_ids = hit_data_ids;

%% Save

save('./fabrizio_longo_2010.mat','fabrizio_longo_2010');

%% Print out

fid = fopen('./fabrizio_longo_2010.txt','w');
write_matrix_file(fid, fabrizio_longo_2010.orfs, fabrizio_longo_2010.ph, fabrizio_longo_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(fabrizio_longo_2010)
end

end

