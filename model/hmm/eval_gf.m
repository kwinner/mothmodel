function ll = eval_gf(f, S)

ll = arrayfun(@(s) sum(f .* (s .^ (0:length(f)-1))), S);

end