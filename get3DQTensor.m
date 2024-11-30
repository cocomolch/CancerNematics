
function [scalarNematicOrder_3D, Q_smooth, n] = get3DQTensor(dECM, mum_per_px, molecularAverageScale, QtensorAverageScale, backgroundIntensity)
px_per_mum = mum_per_px.^(-1);

%first make small gauss filter to correct for discrete pixel errors in orientations
dECM_filt = imgaussfilt3(dECM, 1);
%calc gradients
[dx, dy, dz] = imgradientxyz(dECM_filt);

Ixx = dx.*dx;
Iyy = dy.*dy;
Izz = dz.*dz;
Ixy = dx.*dy;
Iyx = Ixy;
Ixz = dx.*dz;
Izx = Ixz;
Iyz = dy.*dz;
Izy = Iyz;
%% integrate gaussian weighted neighborhood
sigma = ceil(molecularAverageScale.*mum_per_px.^(-1));

Jxx=imgaussfilt3(Ixx,sigma);
Jyy=imgaussfilt3(Iyy,sigma);
Jzz=imgaussfilt3(Izz,sigma);
Jxy=imgaussfilt3(Ixy,sigma);
Jyx = Jxy;
Jxz=imgaussfilt3(Ixz,sigma);
Jzx = Jxz;
Jyz=imgaussfilt3(Iyz,sigma);
Jzy = Jyz;

structureTensor = cell(size(Jxx));
for x = 1:size(Jxx, 1)
    for y = 1:size(Jxx, 2)
        for z = 1:size(Jxx, 3)
            structureTensor{x, y, z} = [Jxx(x, y, z), Jxy(x, y, z), Jxz(x, y, z); Jyx(x,y,z), Jyy(x,y,z),Jyz(x,y,z);Jzx(x,y,z), Jzy(x,y,z), Jzz(x,y,z)];
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

%small coherency and energy regions
toIgnore = orientationEnergy < 0.01 | orientationCoherency < 0.1;


% get Q tensor
nx = cellfun(@(x) x(1), molecularOrientations);
ny = cellfun(@(x) x(2), molecularOrientations);
nz = cellfun(@(x) x(3), molecularOrientations);

Qxx = nx.^2-1/3;
Qyy = ny.^2-1/3;
Qzz = nz.^2-1/3;
Qxy = nx.*ny;
Qxz = nx.*nz;
Qyz = ny.*nz;


Qxx_smooth = imgaussfilt3(Qxx,QtensorAverageScale*px_per_mum );
Qyy_smooth = imgaussfilt3(Qyy,QtensorAverageScale*px_per_mum );
Qzz_smooth = imgaussfilt3(Qzz,QtensorAverageScale*px_per_mum );
Qxy_smooth = imgaussfilt3(Qxy,QtensorAverageScale*px_per_mum );
Qxz_smooth = imgaussfilt3(Qxz,QtensorAverageScale*px_per_mum );
Qyz_smooth = imgaussfilt3(Qyz,QtensorAverageScale*px_per_mum );

%nematic order parameter
Q_smooth = cell(size(Qxx_smooth));
for x = 1:size(Qxx_smooth, 1)
    for y = 1:size(Qxx_smooth, 2)
        for z = 1:size(Qxx_smooth, 3)
            Q_smooth{x, y, z} = [Qxx_smooth(x, y, z), Qxy_smooth(x, y, z), Qxz_smooth(x, y, z); Qxy_smooth(x,y,z), Qyy_smooth(x,y,z),Qyz_smooth(x,y,z);Qxz_smooth(x,y,z), Qyz_smooth(x,y,z), Qzz_smooth(x,y,z)];
        end
    end
end

%"slow syntax" for illustration purposes
%director
n = cellfun(@return3DDirector, Q_smooth, 'UniformOutput', false);
%scalar nematic order
scalarNematicOrder_3D = cellfun(@return3DScalarOrderParameter, Q_smooth);

scalarNematicOrder_3D(toIgnore | isnan(scalarNematicOrder_3D) | dECM_filt < backgroundIntensity) = 0;


