%% Zhang~Schneider, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zhang_schneider_2002.pmid = 12477387;

phenotypes = {'cell size (mean)'; 'cell size (median)'; 'cell size (mode)'};
treatments = {'GYPD'};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Zhang et al supplemental data.xlsx', 'Size Data');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,5:7); % if the dataset is discrete or binary
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% A few corrections
hit_strains(strcmp('TAL004W', hit_strains)) = {'YAL004W'};
hit_strains(strcmp('YELOO1C', hit_strains)) = {'YEL001C'};
hit_strains(strcmp('KL187C', hit_strains)) = {'YKL187C'};
hit_strains(strcmp('YMR41W', hit_strains)) = {'YMR241W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If not possible, eliminate the entry
hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Prepare final dataset

% If the dataset is quantitative:
zhang_schneider_2002.orfs = hit_strains;
zhang_schneider_2002.ph = strcat(phenotypes, '; ', treatments);
zhang_schneider_2002.data = hit_data;

%% Save

save('./zhang_schneider_2002.mat','zhang_schneider_2002');

%% Print out

fid = fopen('./zhang_schneider_2002.txt','w');
write_matrix_file(fid, zhang_schneider_2002.orfs, zhang_schneider_2002.ph, zhang_schneider_2002.data);
fclose(fid);

end

