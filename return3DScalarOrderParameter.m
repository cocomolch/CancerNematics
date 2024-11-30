function S = return3DScalarOrderParameter(Q_smoothEntry)

[~, D] = eig(Q_smoothEntry);
eigenvalues = diag(D);
eigenvaluesSign = sign(eigenvalues);

if sum(eigenvaluesSign) == -1 %two are negative one is positive
    S = 3/2*eigenvalues(eigenvaluesSign==1);
    %alternatively use this
    % S = -3*mean(eigenvalues(eigenvaluesSign==-1));
else
    S = NaN;
end