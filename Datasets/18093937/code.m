%% Szymanski~Goodman, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
szymanski_goodman_2007.pmid = 18093937;

treatments = {'log phase'; 'stationary phase'};

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(szymanski_goodman_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table2.xlsx', 'Sheet1');
data(1,:) = []; % remove column headers

% Get the list of ORFs and the correponding data 
hit_strains = data(:,1);

% Transform  categorical phenotypes into multiple binary phenotypes
[FILENAMES{end+1}, data2] = read_data('xlsread', './raw_data/phenotype_mapping.xlsx','Sheet1');
ph_mapping.ph_orig = data2(2:end,1);
ph_mapping.mat = cell2mat(data2(2:end,2:end));
ph_mapping.ph_new = data2(1,2:end)';

for k = 1 : length(treatments)
    hit_data{k} = zeros(length(hit_strains),length(ph_mapping.ph_new));
    for i = 1 : length(hit_strains)
        ind = find(strcmp(hit_strains{i}, data(:,1)));
        phs = regexp(data{ind,k+1},',','split');
        phs = strtrim(lower(phs));
        for j = 1 : length(phs)
            ind2 = find(strcmp(phs{j}, ph_mapping.ph_orig));
            ind3 = find(abs(ph_mapping.mat(ind2,:))>0);
            hit_data{k}(ind,ind3) = ph_mapping.mat(ind2,ind3);
        end
    end
end

 % Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

phenotypes = strcat({'lipid droplets '}, ph_mapping.ph_new);
conditions = [repmat({'log phase'}, length(phenotypes),1); repmat({'stationary phase'}, length(phenotypes),1)];
phenotypes = [phenotypes; phenotypes];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [170 701 702 703 704 705 706 707 5387 5388 5389 5390 5391 5392 5393 5394]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
szymanski_goodman_2007.orfs = hit_strains;
szymanski_goodman_2007.ph = hit_data_names;
szymanski_goodman_2007.data = [hit_data{1} hit_data{2}];
szymanski_goodman_2007.dataset_ids = hit_data_ids;

% Remove the phenotypes/datasets that (as a result of the conversion 
% from categorical to binary phenotypes) ended up not having any hits

ix = find(sum(abs(szymanski_goodman_2007.data),1)==0);
szymanski_goodman_2007.ph(ix) = [];
szymanski_goodman_2007.data(:,ix) = [];
szymanski_goodman_2007.dataset_ids(ix) = [];

%% Save

save('./szymanski_goodman_2007.mat','szymanski_goodman_2007');

%% Print out

fid = fopen('./szymanski_goodman_2007.txt','w');
write_matrix_file(fid, szymanski_goodman_2007.orfs, szymanski_goodman_2007.ph, szymanski_goodman_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(szymanski_goodman_2007)
end

end

