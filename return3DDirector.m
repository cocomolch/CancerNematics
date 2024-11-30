function n = return3DDirector(Q_smoothEntry)

 [V, D] = eig(Q_smoothEntry);
eigenvalues = diag(D);
eigenvaluesSign = sign(eigenvalues);

if sum(eigenvaluesSign) == -1 %two are negative one is positive
    S = 3/2*eigenvalues(eigenvaluesSign==1);
    n = V(:, eigenvaluesSign==1);
    %alternatively use this
    % S = -3*mean(eigenvalues(eigenvaluesSign==-1));
else
    n = NaN;
    S = NaN;
end