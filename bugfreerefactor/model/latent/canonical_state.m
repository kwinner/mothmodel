function Q = canonical_state( y, N )

K = numel(y);

%running total
nAlive = 0;

%birth and death records
births = zeros(1, K+1);
deaths = zeros(1, K+1);

%the actual result
Q      = zeros(K+1, K+1);

for i = 1:K
	if y(i) >= nAlive 
		%need to birth individuals to explain y(i)
		delta = y(i) - nAlive;
		births(i) = births(i) + delta;
		nAlive = nAlive + delta;
	else
		%too many individuals alive, kill some off
		quota = nAlive - y(i);
		for j = 1:i
			if deaths(j) < births(j)
				delta = min(births(j) - deaths(j), quota); %how many to kill off
				deaths(j) = deaths(j) + delta;
				nAlive = nAlive - delta;
				quota  = quota - delta;
				Q(j,i) = Q(j,i) + delta;

				%stop this senseless killing
				if quota <= 0; break; end
			end
		end
	end
end

%kill any remaining survivors
for i = 1:K+1
	if deaths(i) < births(i)
		delta = births(i) - deaths(i);
		Q(i,K+1) = Q(i,K+1) + delta;
	end
end

%balance all unobserved individuals along the diagonal
unobserved = N - sum(Q(:));
even   = fix(unobserved / (K+1));
uneven = rem(unobserved,  K+1);

if even ~= 0
	Q = Q + even .* eye(K+1,K+1);
end
if uneven ~= 0
	Q(:, 1:uneven) = Q(:, 1:uneven) + eye(K+1, uneven);
end

end

