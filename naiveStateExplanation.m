function [ q, n ] = naiveStateExplanation( y, N )
%NAIVESTATEEXPLANATION Return the "simplest" q (and corresponding n) explaining a set of observations y

T = numel(y);

nAlive = 0;
nDead  = 0;
births = zeros(1,T+1);
deaths = zeros(1,T+1);
q      = zeros(T+1,T+1);

for i = 1:T
    if y(i) >= nAlive %more individuals needed to explain y(i)
        %birth(I(i-1),y(i)-#alive)
        babies = y(i) - nAlive;
        births(i) = births(i) + babies;
        nAlive = nAlive + babies;
    else %too many individuals alive, kill some
        %kill(I(i-1),y(i)-nAlive)
        quota = nAlive - y(i);
        for j = 1:i
            if deaths(j) < births(j)
                avail = births(j) - deaths(j); %how many we can kill from j
                victims = min(quota,avail);    %how many we actually need to kill
                %kill them:
                deaths(j) = deaths(j) + victims;
                nAlive    = nAlive - victims;
                nDead     = nDead + victims;
                quota     = quota - victims;
                q(j,i)    = q(j,i) + victims;
                %do we need to keep killing??
                if quota <= 0; break; end
            end
        end
    end
end
nObserved = nAlive + nDead;
%kill all the final survivors
for i = 1:T+1
    if deaths(i) < births(i)
        victims = births(i) - deaths(i);
        q(i,T+1) = q(i,T+1) + victims;
        nAlive = nAlive - victims;
        nDead  = nDead + victims;
        if nAlive <= 0; break; end
    end
end
%stick all the 'unobserved' individuals in at the beginning: q(1,1)
q(1,1) = N - sum(q(:));

n = abundancy(q);

end

