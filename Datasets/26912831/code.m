%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
saka_kobayashi_2016.pmid = 26912831;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(saka_kobayashi_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/geldata_201604042348.csv', 'Delimiter',',');

% Get the list of ORFs and the correponding data 
hit_strains = data.ORFName;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Get the data itself
hit_data = nan(size(hit_strains,1),2);

% Adjust the data
inds = find(strcmp('N.D.', data.rDNAStability));
data.rDNAStability(inds) = {'NaN'};
data.rDNAStability = cellfun(@str2num, data.rDNAStability);
hit_data(:,1) = -data.rDNAStability+2;  % wt = 0, h = negative values

inds = find(strcmp('N.D.', data.rDNACopyNumber) | strcmp('-', data.rDNACopyNumber));
data.rDNACopyNumber(inds) = {'NaN'};
data.rDNACopyNumber = strtrim(data.rDNACopyNumber);

% Average the multiple calls
inds = find(~cellfun(@isempty, regexp(data.rDNACopyNumber, ',')));
for i = 1 : length(inds)
    tmp = regexp(data.rDNACopyNumber{inds(i)}, ',', 'split');
    hit_data(inds(i),2) = mean(cellfun(@str2num, tmp));
end
inds = find(cellfun(@isempty, regexp(data.rDNACopyNumber, ',')));
hit_data(inds,2) = cellfun(@str2num, data.rDNACopyNumber(inds));
hit_data(:,2) = hit_data(:,2)-2;  % wt = 0, s = positive values

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [756 746];
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
saka_kobayashi_2016.orfs = hit_strains;
saka_kobayashi_2016.ph = hit_data_names;
saka_kobayashi_2016.data = hit_data;
saka_kobayashi_2016.dataset_ids = hit_data_ids;

%% Save

save('./saka_kobayashi_2016.mat','saka_kobayashi_2016');

%% Print out

fid = fopen('./saka_kobayashi_2016.txt','w');
write_matrix_file(fid, saka_kobayashi_2016.orfs, saka_kobayashi_2016.ph, saka_kobayashi_2016.data);
fclose(fid);

end

