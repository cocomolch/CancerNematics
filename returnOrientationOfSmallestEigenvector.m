function orientation = returnOrientationOfSmallestEigenvector(structureTensorEntry, casy)

[V, D] = eig(structureTensorEntry);
eigenvalues = diag(D);
[~, sortInd] = sort(eigenvalues, 'ascend');

smallestEigenvector = V(:, sortInd(1));
orientationEnergy = abs(eigenvalues(sortInd(end))) + abs(eigenvalues(sortInd(1)));
if length(eigenvalues) == 3
    orientationCoherency = abs(eigenvalues(sortInd(3)) - eigenvalues(sortInd(1)))./abs(eigenvalues(sortInd(3))) + abs(eigenvalues(sortInd(1)));
elseif length(eigenvalues) == 2
    orientationCoherency = abs(eigenvalues(1) - eigenvalues(2))./abs(eigenvalues(1)) + abs(eigenvalues(2));
else
    print('DIMENSION ERROR')
end
switch casy%very slow but for illustration purposes
    case 1
        orientation = smallestEigenvector;
    case 2
        orientation = orientationCoherency;
    case 3
        orientation = orientationEnergy;
end