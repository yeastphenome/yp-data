%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
saleem_aitchison_2010.pmid = 20694151;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(saleem_aitchison_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pone.0011953.s004.xlsx', '1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(6:end,2);

% Get the data itself
hit_data = data(6:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

hit_strains(find(strcmp('YYKL138C', hit_strains))) = {'YKL138C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data2 = zeros(length(hit_strains), 5);
hit_data(~cellfun(@isnumeric, hit_data(:,1)),1) = {NaN};
hit_data2(:,1) = cell2mat(hit_data(:,1));
for i = 1 : length(hit_strains)
    switch hit_data{i,2}
        case 'A'
            % nothing to do, phenotype is normal
        case 'B'
            hit_data2(i,2) = -1;
        case 'C'
            hit_data2(i,2) = 1;
        case 'D'
            hit_data2(i,3) = -1;
        case 'E'
            hit_data2(i,3) = 1;
        case 'G'
            % not clear about what phenotype this is, skipping
        case 'H'
            % not clear about what phenotype this is, skipping
        case 'J'
            hit_data2(i,4) = -1;
        case 'K'
            hit_data2(i,5) = -1;
        case 'L'
            hit_data2(i,5) = -0.5;
    end
    
end

hit_data = hit_data2;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11830 11831 11832 11834 11835]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
saleem_aitchison_2010.orfs = hit_strains;
saleem_aitchison_2010.ph = hit_data_names;
saleem_aitchison_2010.data = hit_data;
saleem_aitchison_2010.dataset_ids = hit_data_ids;

%% Save

save('./saleem_aitchison_2010.mat','saleem_aitchison_2010');

%% Print out

fid = fopen('./saleem_aitchison_2010.txt','w');
write_matrix_file(fid, saleem_aitchison_2010.orfs, saleem_aitchison_2010.ph, saleem_aitchison_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(saleem_aitchison_2010)
end

end

