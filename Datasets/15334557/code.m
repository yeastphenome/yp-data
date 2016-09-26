%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mollapour_piper_2004.pmid = 15334557;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mollapour_piper_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('textscan','./raw_data/Table1.txt', '%s');
[FILENAMES{end+1}, data2] = read_data('textscan','./raw_data/Table2.txt', '%s');
[FILENAMES{end+1}, data3] = read_data('textscan','./raw_data/Table3.txt', '%s');

% Get the list of ORFs and the correponding data 
hit_strains1 = clean_orf(data1);
hit_strains2 = clean_orf(data2);
hit_strains3 = clean_orf(data3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds)); 

hit_strains2(ismember(hit_strains2, {'TMR205C'})) = {'YMR205C'};
hit_strains2(ismember(hit_strains2, {'YMLO48W'})) = {'YML048W'};
hit_strains2(ismember(hit_strains2, {'YLRO25W'})) = {'YLR025W'};
hit_strains2(ismember(hit_strains2, {'YKLO54W'})) = {'YKL054W'};
hit_strains2(ismember(hit_strains2, {'YCLO58C'})) = {'YCL058C'};

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds)); 

inds = find(~is_orf(hit_strains3));
disp(hit_strains3(inds)); 

hit_strains = unique([hit_strains1; hit_strains2; hit_strains3]);
hit_data = zeros(length(hit_strains),3);

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1) = -1;
[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,2) = -1;
[~,ind1,ind2] = intersect(hit_strains, hit_strains3);
hit_data(ind1,3) = 1;


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [487; 554; 555];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
mollapour_piper_2004.orfs = hit_strains;
mollapour_piper_2004.ph = hit_data_names;
mollapour_piper_2004.data = hit_data;
mollapour_piper_2004.dataset_ids = hit_data_ids;

%% Save

save('./mollapour_piper_2004.mat','mollapour_piper_2004');

%% Print out

fid = fopen('./mollapour_piper_2004.txt','w');
write_matrix_file(fid, mollapour_piper_2004.orfs, mollapour_piper_2004.ph, mollapour_piper_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mollapour_piper_2004)
end

end

