
target_folder = 'first_experiment';
base_folder = '/Users/gbernstein/Documents/Work/MLDS/Moth/coverage_results';
folder = [base_folder, filesep, target_folder, filesep];
files = dir([folder, '*.out']);

%lines in lines var
mu_line = 2;
sigma_line = 3;
lambda_line = 4;
T_line = 5;
alpha_line = 6;
N_line = 7;
n_iters_line = 8;

alg_strs = {'zonn','gp'};
starts =  [10, 16];
offsets = [0,1,2,3];
mean_strs = {'theta','CI','error','coverage'};
means = zeros(length(starts),length(offsets),4);
param_strs = {'mu','sigma','lambda','T','u','alpha','N','n_iters'};
target_strs = {'mu','sigma','lambda','N'};

%fill out empty lines
to_print = cell(length(files)+2, length(param_strs) + length(target_strs)*length(alg_strs)*length(mean_strs));
to_print(2,1:length(param_strs)) = param_strs;
alg_inds = repmat([1,2],1,4*length(mean_strs));
mean_inds = repmat([1,1,2,2,3,3,4,4],1,4);
param_inds = repmat(1:length(param_strs),length(mean_strs)*length(alg_strs));
to_print(2,length(param_strs)+1:end) = arrayfun(@(i) [alg_strs{alg_inds(i)},'_',mean_strs{mean_inds(i)},'_',param_strs{param_inds(i)}],1:length(mean_inds),'UniformOutput',false);
params = struct();
for f = 1:length(files)
    
    %split file into lines
    str = fileread([folder,filesep,files(f).name]);
    lines = strsplit(str,'\n')';
    
    %mu
    temp = strsplit(lines{mu_line},' ');
    params.mu = str2double(temp{2});    
    %sigma
    temp = strsplit(lines{sigma_line},' ');
    params.sigma = str2double(temp{2});    
    %lambda
    temp = strsplit(lines{lambda_line},' ');
    params.lambda = str2double(temp{2});    
    %T and u
    temp = strsplit(lines{T_line},' ');
    temp(1) = [];
    temp{1} = temp{1}(2:end);
    temp{end} = temp{end}(1:end-1);
    params.T = cellfun(@str2double,temp);    
    params.u = mean(diff(params.T));
    params.T =  strjoin(arrayfun(@num2str,params.T,'UniformOutput',false),';');
    %alpha
    temp = strsplit(lines{alpha_line},' ');
    params.alpha = str2double(temp{2});    
    %N
    temp = strsplit(lines{N_line},' ');
    params.N = str2double(temp{2});    
    %n_iters
    temp = strsplit(lines{n_iters_line},' ');
    params.n_iters = str2double(temp{2});

    %read in means
    for s = 1:length(starts)
        for o = 1:length(offsets)
            temp = strsplit(lines{starts(s) + offsets(o)},'=');
            temp = strsplit(temp{2},' ');
            temp(1) = [];
            temp{1} = temp{1}(2:end);
            temp{end} = temp{end}(1:end-1);
            means(s,o,:) = cellfun(@str2double,temp);
        end
    end
    
    %fill out to_print
    i = 1;
    for p = 1:length(param_strs)
        to_print{f+2,i} = params.(param_strs{p});
        i = i + 1;
    end
    for t = 1:length(target_strs)
        for o = 1:length(mean_strs)
            for s = 1:length(alg_strs)            
                to_print{f+2,i} = means(s,o,t);
                i = i + 1;
            end
        end
    end
end

%print to file
fid = fopen([folder, filesep, 'aaa_table.csv'],'w+');
for r = 1:size(to_print,1)
    for c = 1:size(to_print,2)
        if r <= 2
            fprintf(fid,'%s,',to_print{r,c});
        else
            if c == find(strcmp(param_strs,'T'))  
                fprintf(fid,'%s,',to_print{r,c});
            else
                fprintf(fid,'%.2f,',to_print{r,c});
            end
            
        end
    end
    fprintf(fid,'\n');
end

    
    
    
    
    
