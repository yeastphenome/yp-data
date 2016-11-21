%% Junne~Hoepfner, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
junne_hoepfner_2015.pmid = 25616894;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(junne_hoepfner_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

spreadsheets = {'Decatransin Cpd 1 Exp 1','Decatransin Cpd 2 Exp 2','Cotransin Cpd 2','Cotansin Cpd 3 (HUN-7293)'};

for d = 1 : length(spreadsheets)

    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Junnet_et_al_HIPHOPrawdata.xlsx', spreadsheets{d});

    % Get the list of ORFs and the correponding data 

    hom_indx = find(strcmp('HOP',data(:,1)));
    het_indx = find(strcmp('HIP',data(:,1)));
    hom_strains{d} = data(hom_indx,5);
    het_strains{d} = data(het_indx,5);

    hom_data{d} = cell2mat(data(hom_indx, 7));
    het_data{d} = cell2mat(data(het_indx, 7));

    % Eliminate all white spaces & capitalize
    hom_strains{d} = clean_orf(hom_strains{d});
    het_strains{d} = clean_orf(het_strains{d});

    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hom_strains{d}));
    hom_strains{d}(inds) = [];
    hom_data{d}(inds) = [];

    inds = find(~is_orf(het_strains{d}));
    het_strains{d}(inds) = [];
    het_data{d}(inds) = [];

end

%% Combine the data into one matrix

all_strains = unique([vertcat(hom_strains{:}); vertcat(het_strains{:})]);
all_data = nan(length(all_strains), 8);

for i = 1 : 4
    
    [~, ind1, ind2] = intersect(all_strains, hom_strains{i});
    all_data(ind1, i*2-1) = hom_data{i}(ind2);

    [~, ind1, ind2] = intersect(all_strains, het_strains{i});
    all_data(ind1, i*2) = het_data{i}(ind2);

end

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [772; 773; 5262; 5263; 5264; 5265; 5266; 5267];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
junne_hoepfner_2015.orfs = all_strains;
junne_hoepfner_2015.ph = hit_data_names;
junne_hoepfner_2015.data = all_data;
junne_hoepfner_2015.dataset_ids = hit_data_ids;

%% Save

save('./junne_hoepfner_2015.mat','junne_hoepfner_2015');

%% Print out

fid = fopen('./junne_hoepfner_2015.txt','w');
write_matrix_file(fid, junne_hoepfner_2015.orfs, junne_hoepfner_2015.ph, junne_hoepfner_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(junne_hoepfner_2015)
end

end

