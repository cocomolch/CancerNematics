function [S, nx, ny] = get2DQtensor(dECM_2D, mum_per_px_2D, molecularAverageScale, QtensorAverageScale, filterOrientations,orientationEnergyThreshold,orientationCoherencyThreshold, backgroundIntensity)


%first make small gauss filter to correct for pixel errors in orientations
dECM_filt = imgaussfilt(dECM_2D, 1);
%calc gradients
[dx, dy] = imgradientxy(dECM_filt);

Ixx = dx.*dx;
Iyy = dy.*dy;
Ixy = dx.*dy;
Iyx = Ixy;
%% integrate gaussian weighted neighborhood
sigma = floor(molecularAverageScale.*mum_per_px_2D.^(-1));

Jxx=imgaussfilt(Ixx,sigma);
Jyy=imgaussfilt(Iyy,sigma);
Jxy=imgaussfilt3(Ixy,sigma);
Jyx = Jxy;


structureTensor = cell(size(Jxx));
for x = 1:size(Jxx, 1)
    for y = 1:size(Jxx, 2)
        for z = 1:size(Jxx, 3)
            structureTensor{x, y, z} = [Jxx(x, y, z), Jxy(x, y, z); Jyx(x,y,z), Jyy(x,y,z)];
        end
    end
end
clear Jxx Jyy Jzz Jxy Jyx Jyz Jzy Jxz Jzx Ixx Iyy Izz Ixy Iyx Iyz Izy Ixz Izx

% now evaluate eigensystem at each pixel ("slow syntax" for illustration purposes)

molecularOrientations = cellfun(@(x)returnOrientationOfSmallestEigenvector(x, 1), structureTensor, 'UniformOutput',false);
orientationCoherency = cellfun(@(x)returnOrientationOfSmallestEigenvector(x, 2), structureTensor, 'UniformOutput',false);
orientationCoherency = cell2mat(orientationCoherency);
orientationEnergy = cellfun(@(x)returnOrientationOfSmallestEigenvector(x, 3), structureTensor, 'UniformOutput',false);
orientationEnergy = cell2mat(orientationEnergy);


% get Q tensor
nx = cellfun(@(x) x(1), molecularOrientations);
ny = cellfun(@(x) x(2), molecularOrientations);

Qxx=nx.^2-1/2;
Qxy=nx.*ny;

Qxx_smooth=imgaussfilt(Qxx,QtensorAverageScale*mum_per_px_2D^(-1));
Qxy_smooth=imgaussfilt(Qxy,QtensorAverageScale*mum_per_px_2D^(-1));
S=2*sqrt((Qxx_smooth.^2+Qxy_smooth.^2)); %nematic order parameter

if filterOrientations
    %small coherency and energy regions
    toIgnore = orientationEnergy < orientationEnergyThreshold | orientationCoherency < orientationCoherencyThreshold;
    S(toIgnore | dECM_filt < backgroundIntensity) = NaN;
end

