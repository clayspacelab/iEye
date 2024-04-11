function i = moditer(i, nloop)
% Allows use of mod for iteration
    i = i - 0.1;
    i = mod(i,nloop);
    i = round(i);
end