function start_experiment(name)
%prepare a results directory and a meta file for this experiment

%form a timestamp
timestamp = now;

%create the output directory
resultsDir = fullfile(getenv('RESULTSDIR') ...
	                  ,[datestr(floor(timestamp)) ...
	                    , ' ' ...
	                    , name]);
%while the results directory already exists, update the folder name to avoid collision
while exist(resultsDir, 'file') == 7
	
end

end

