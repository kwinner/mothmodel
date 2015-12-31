%set default settings for moth model experiments

% T     = [0,8,20];  %observation times
% N_hat = 20;    %mean superpopulation size
% n_max = N_hat; %cap on abundance at each time

% alpha   = 1.0;   %detection probability
% mu      = 8; %mean arrival
% sigma   = 4; %var of arrival
% lambda  = 3; %mean service

% arrival rate function
% arrivalDistn = makedist('Normal', 'mu', mu, 'sigma', sigma);
% rateFunc     = @arrivalDistn.pdf;

% %service distribution
% serviceDistn = makedist('Exp', 'mu', lambda);


% [S,Z] = sample_pop(arrivalDistn, serviceDistn, N_hat);
% n = abundance2(S,Z,T);
% y = sample_obs(alpha,n);

% P = equivalence_distn(arrivalDistn, serviceDistn, T);
% Q_0 = hist3([S,S+Z], 'Edges', {[T inf], [T inf]});

% [y,n,S,Z,Q_0] = samplestate(mu,sigma,lambda,N_hat,alpha,T);
% P = birthdeath_pmf(mu,sigma,lambda,N_hat,T);

% [~,Q_history] = mcmc(y, Q_0, P, N_hat, alpha, T, 'nIterations', 100000);

% Q = Q_history{end};

% how many of each q bin we've filled so far
% T_plot = -5:.1:45;
% abundances = zeros(10000, length(T_plot));
% for i = 1:10000
% 	if mod(i-1, 10) == 0
% 		disp(i)
% 	end
% 	counts = zeros(size(Q_history{1000+9*(i-1)}));
% 	Q = Q_history{1000+9*(i-1)};

% 	% list of birth and death times
% 	draws = zeros(N_hat,2); % [birth time, death time]
% 	iDraw = 1;

% 	%append -inf and inf to t
% 	tWithBounds = [-inf, T, inf];

% 	for iBirth = 1:size(Q,1)
%     	%create a truncated birth distribution which forces births to be in the right window
%     	smin = tWithBounds(iBirth);
%     	smax = tWithBounds(iBirth + 1);
% 	    birthDist = truncate(makedist('Normal', 'mu', mu, 'sigma', sigma), smin, smax);

% 	    for iDeath = iBirth:size(Q,2)
% 	        %truncate the lifespan distribution pessimistically
% 	        %min lifespan = mindeathtime - maxbirthtime
% 	        zmin = tWithBounds(iDeath) - smax;
% 	        %max lifespan = maxdeathtime - minbirthtime
% 	        zmax = tWithBounds(iDeath + 1) - smin;

% 	        lifespanDist = truncate(makedist('Exponential', 'mu', lambda), zmin, zmax);

% 	        while counts(iBirth, iDeath) < Q(iBirth, iDeath)
% 	        	% disp(sprintf('[%d,%d]: %d/%d', iBirth, iDeath, counts(iBirth, iDeath), Q(iBirth, iDeath)))
% 	            %draw a birth and death
% 	            birth    = random(birthDist);
% 	            lifespan = random(lifespanDist);
% 	            death    = birth + lifespan;

% 	            %check if the individual is needed in the counts
% 	            myDeath = sum(T < death) + 1;
% 	            if counts(iBirth, myDeath) < Q(iBirth, myDeath)
% 	                draws(iDraw, :) = [birth, lifespan];
% 	                counts(iBirth, myDeath) = counts(iBirth, myDeath) + 1;
% 	                iDraw = iDraw + 1;
% 	            end
%             end
%         end
%     end

%     abundances(i,:) = abundance2(draws(:,1), draws(:,2), T_plot);
% end

% bincounts = histc(abundances, 0:1:max(abundances(:)));

%plot w/ conf intervals
figure('units','normalized','position',[.1 .1 .8 .8]); hold on
% plot([0,45],[0,0],'Color','black')
back = fill([T_plot'; flipdim(T_plot', 1)], [mean(abundances)' + 2.*std(abundances)'; flipdim(mean(abundances)' - 2.*std(abundances)', 1)], [7 7 7]/8);
back.FaceAlpha = 0.8;
back.EdgeAlpha = 0.8;
back.EdgeColor = [7 7 7]/12;
plot(T_plot, mean(abundances), 'LineWidth', 2)
title('Abundance over time', 'FontSize', 16)
xlabel('Days after May 1st', 'FontSize', 14)
ylabel('Abundance', 'FontSize', 14)
set(gca, 'FontSize', 12)

figure('units','normalized','position',[.1 .1 .8 .8]);
imagesc(flipud(log(bincounts)))
ax = gca;
colormap(ax, 'default')
ax.XTickLabel = 0:5:45;
ax.YTickLabel = 32:-5:2;
title('Abundance over time', 'FontSize', 16)
xlabel('Days after May 1st', 'FontSize', 14)
ylabel('Abundance', 'FontSize', 14)
set(ax, 'FontSize', 12)

figure('units','normalized','position',[.1 .1 .8 .8]);
imagesc(flipud(log(bincounts)))
ax = gca;
colormap(gca, 'hot')
ax.XTickLabel = 0:5:45;
ax.YTickLabel = 32:-5:2;
title('Abundance over time', 'FontSize', 16)
xlabel('Days after May 1st', 'FontSize', 14)
ylabel('Abundance', 'FontSize', 14)
set(ax, 'FontSize', 12)

figure('units','normalized','position',[.1 .1 .8 .8]);
imagesc(flipud(log(bincounts)))
ax = gca;
colormap(gca, 'jet')
ax.XTickLabel = 0:5:45;
ax.YTickLabel = 32:-5:2;
title('Abundance over time', 'FontSize', 16)
xlabel('Days after May 1st', 'FontSize', 14)
ylabel('Abundance', 'FontSize', 14)
set(ax, 'FontSize', 12)
