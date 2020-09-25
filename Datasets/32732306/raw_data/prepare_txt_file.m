s3 = load('raw_data/File_S3.mat');

colony_fitness = s3.sumint_40.prim.norm_fitness_med;
lagVstall = s3.sumint_40.prim.tdev_med;

orfs = s3.sumint_40.uorfs;

fid = fopen('raw_data/File_S3_data.txt','w');
for i = 1 : length(orfs)
    fprintf(fid, '%s\t%.3f\t%.3f\t%.3f\t%.3f\n', orfs{i}, colony_fitness(i,1), colony_fitness(i,2), lagVstall(i,1), lagVstall(i,2));
end
fclose(fid);