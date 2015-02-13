function [zmin, zmax] = lifespan_domain (birth, death_min, death_max)
%note: birth can be a vector, death_min/max should be scalars

%zmin is the minimum lifetime to live from birth to death_min
zmin = death_min - birth;
zmin(death_min == birth) = 0; %this is to cover when death_min, birth == -inf
zmin = max(zmin, 0);          %lifetime cannot be negative, so the true min is 0

%zmax is the maximum lifetime to live from birth to death_max
zmax = death_max - birth;
zmax(death_max == birth) = inf; %this is to cover when death_max, birth == inf
zmax = max(zmax, zmin);         %keep zmax >= zmin, usually only if zmin brought up to 0

end