%% Hoon~Nislow, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hoon_nislow_2011.pmid = 22384333;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hoon_nislow_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/TableS1.xlsx');
hit_strains1 = data1(5:end,2);
hit_data1 = cell2mat(data1(5:end,4));

hit_strains1 = clean_orf(hit_strains1);
hit_strains1 = translate(hit_strains1);

inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});


% Load file
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/TableS2.xlsx');
hit_strains2 = data2(5:end,2);
hit_data2 = cell2mat(data2(5:end,4));

hit_strains2 = clean_orf(hit_strains2);
hit_strains2 = translate(hit_strains2);

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});


% Load file
[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/TableS3.xlsx');

% Get the list of ORFs
hit_strains3 = data3(5:end, 1);
hit_data3 = -cell2mat(data3(5:end, 6)); 

% Clean up ORFs
hit_strains3 = cellfun(@(x) strtok(x, ':'), hit_strains3, 'UniformOutput', false); 
hit_strains3 = clean_orf(hit_strains3);

hit_strains3 = translate(hit_strains3);

inds = find(~is_orf(hit_strains3));
disp(hit_strains3(inds)); 

[hit_strains3, hit_data3] = grpstats(hit_data3, hit_strains3, {'gname','mean'});


hit_strains = unique([hit_strains1; hit_strains2; hit_strains3]);
hit_data = zeros(length(hit_strains),3);

[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1) = hit_data1(ind1);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,2) = hit_data2(ind1);
[~,ind1,ind2] = intersect(hit_strains3, hit_strains);
hit_data(ind2,3) = hit_data3(ind1);


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [519; 755; 16253];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hoon_nislow_2011.orfs = hit_strains;
hoon_nislow_2011.ph = hit_data_names;
hoon_nislow_2011.data = hit_data;
hoon_nislow_2011.dataset_ids = hit_data_ids;

%% Save

save('./hoon_nislow_2011.mat','hoon_nislow_2011');

%% Print out

fid = fopen('./hoon_nislow_2011.txt','w');
write_matrix_file(fid, hoon_nislow_2011.orfs, hoon_nislow_2011.ph, hoon_nislow_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hoon_nislow_2011)
end

end