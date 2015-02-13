function [ data ] = mpg ( directory )

NUM_TOP_SPECIES = 5;

%read data
if ~exist('directory', 'var') || exist(directory, 'file') ~= 7
	directory = '/Users/kwinn/Work/Data/moth/';
end

filename_pollard = fullfile(directory, 'MPGPollard.csv');
filename_timed   = fullfile(directory, 'MPGTimed.csv');

data_pollard = csv2struct(filename_pollard, '%d%d%s%s%s%d%d%d%d%d%d');
data_timed   = csv2struct(filename_timed,   '%d%d%s%s%s%s%d%d%d%d%d');

%pick out the top K species by total observation count
species = union(unique({data_pollard.Species}), unique({data_timed.Species}));
species = species(cellfun(@(x) isempty(x), strfind(species,'Unknown'), 'UniformOutput', true));
species = species(cellfun(@(x) isempty(x), strfind(species,'Common/Western'), 'UniformOutput', true));

species_totals = zeros(size(species));
for ispecies = 1:numel(species)
	species_totals(ispecies) = species_totals(ispecies) + ...
							   sum([data_pollard(strcmp({data_pollard.Species}, species{ispecies})).Total]) + ...
							   sum([data_timed(  strcmp({data_timed.Species},   species{ispecies})).Max]);
end

[species_totals, isort] = sort(species_totals, 'descend');
species = species(isort);

% bar(species_totals);
% xticklabel_rotate(1:numel(species), 90, species);

species        = species(1:NUM_TOP_SPECIES);
species_totals = species_totals(1:NUM_TOP_SPECIES);

%convert the date in all data to a delta from first observation
date_pollard = datenum({data_pollard.survey_date}, 'mm/dd/yy');
date_timed   = datenum({data_timed.survey_date},   'mm/dd/yy');

t_0 = min([date_pollard; date_timed]);

date_pollard = date_pollard - t_0;
date_timed   = date_timed   - t_0;

for irecord = 1:numel(data_pollard)
	data_pollard(irecord).survey_date = date_pollard(irecord);
end

for irecord = 1:numel(data_timed)
	data_timed(irecord).survey_date = date_timed(irecord);
end

%group sites by 'lump'
lumps = {}; lumptitles = {};
filename_lumps = fullfile(directory, 'lumps_binary.csv');
lumps_file = fopen(filename_lumps);
lump_raw = fgets(lumps_file);
while(ischar(lump_raw))
	lumpsplit = strtrim(strsplit(lump_raw, ','));
	lumptitles{end+1} = lumpsplit{1};
	lumps{end+1} = lumpsplit(2:end);

	lump_raw = fgets(lumps_file);
end
fclose(lumps_file);

for jspecies = 1:numel(species)
	% figure
	for ilump = 1:numel(lumps)
		obs_ij_pollard = data_pollard(and(strcmp({data_pollard.Species}, species{jspecies}), ismember({data_pollard.Route}, lumps{ilump})));
		obs_ij_timed   = data_timed(  and(strcmp({data_timed.Species  }, species{jspecies}), ismember({data_timed.Route},   lumps{ilump})));

		% y_raw = [obs_ij_timed.Max];
		% t_raw = [obs_ij_timed.survey_date];
		y_raw = [obs_ij_pollard.Total,       obs_ij_timed.Max];
		t_raw = [obs_ij_pollard.survey_date, obs_ij_timed.survey_date];

		%sum up all counts from the same day to get the final obs times and total counts
		t_unique = unique(t_raw); T = numel(t_unique);
		t = zeros(1,T);
		y = zeros(1,T);
		for it = 1:numel(t_unique)
			t(it) = t_unique(it);
			y(it) = sum(y_raw(t_raw == t_unique(it)));
		end

		% subplot(numel(lumps), numel(species), jspecies + (ilump-1)*numel(species)); plot(t,y,'o');
		% subplot(numel(lumps), 1, ilump); plot(t,y,'o');
		% xlim([0,150]);
		% if jspecies == 1
		% 	ylabel(lumptitles{ilump})
		% end
		% if ilump == numel(lumps)
		% 	xlabel(species{jspecies})
		% end

		theta = MLtheta(y, struct('t', t, 'alpha', 1),'objective', 'zonn', 'theta_0', [mean(t), 1/mean(t), 7, mean(y)/.125], 'theta_ub', [max(t), mean(t), max(t)-min(t), sum(y)]);

		t_plot = min(t):.1:max(t);
		p = ppdf(theta, struct('t', t_plot, 'alpha', .25, 'N', theta.N));
		pt = arrayfun(@(k) params.alpha.*sum(sum(p(1:k, k+1:numel(t_plot)+1))), 1:t_plot);

		hold on
		plot(t_plot, pt)

		keyboard
	end
end

end

