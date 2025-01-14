classdef DefectAnalysis
    %ECM_ALIGNMENT Measures the alignment of the ECM
    
    properties
        nematicParameterThreshold;
        defectCharge; %defects with charge +-0.5+-chargeDelta are considered
        chargeDelta;
        radiusLineIntegralRange;%typically just one value (width of kernel path)
        nanRegionSizeToAllow;% forbidden area size over which we just coarse grain
        mergeDefectsLengthScale;%if we want to merge defects for some reason (legacy)
    end
    
    
    methods (Static)
        
        function ECMDefects = findDefects(ECMOrientations, nx,ny,S, radiusLineIntegralRange, nematicParameterThreshold, defectCharge, chargeDelta, nanRegionSizeToAllow, mergeDefectsLengthScale)
            %% ECMOrientations = list of coordinates and associated directions ([X, Y] = meshgrid(1:size(nx, 1), 1:size(ny, 2)); ECMOrientations = [X(:), Y(:), nx_smooth(:), ny_smooth(:)])
            
            ECMDefects.nematicParameterThreshold    = nematicParameterThreshold;
            ECMDefects.defectCharge                 = defectCharge;
            ECMDefects.chargeDelta                  = chargeDelta;
            ECMDefects.radiusLineIntegralRange      = radiusLineIntegralRange;
            ECMDefects.mergeDefectsLengthScale      = mergeDefectsLengthScale;
            %% get gradient map of angles
            thetaMap = acos(nx);
            %           thetaMap = rad2deg(thetaMap);
            [theta_x, theta_y] = gradient(thetaMap);
            %setting delta theta to zeros whgere NaNs are
            theta_x(isnan(theta_x)) = 0;
            theta_y(isnan(theta_y)) = 0;
            %% make NaN-Map in order to delete defetcs in vicinity of undefined areas
            nanMap = DefectAnalysis.makeNanMap(thetaMap, ECMOrientations, nanRegionSizeToAllow);
            
            %init all defect variables
            positiveCentroids_all           = [];
            negativeCentroids_all           = [];
            positiveCharges_all             = [];
            negativeCharges_all             = [];
            positiveAngle_all               = [];
            negativeAngle_all               = [];
            positiveNematicParameter_all    = [];
            negativeNematicParameter_all    = [];
            numberPosDefectsAtLengthScale  = zeros(1, length(radiusLineIntegralRange));
            numberNegDefectsAtLengthScale  = zeros(1, length(radiusLineIntegralRange));
            
            for loopInegralSizeIndex = 1:length(radiusLineIntegralRange)
                currentRadiusLineIntegral = radiusLineIntegralRange(loopInegralSizeIndex);
                %% make ring or square kernel
                 [deltaX, deltaY, ring_kernel]   = DefectAnalysis.makeRingKernel(round(currentRadiusLineIntegral/2));
                 %[deltaX,deltaY,ring_kernel] = DefectAnalysis.makeSquareKernel(round(currentRadiusLineIntegral/2));

                %% define dtheta/dx and dtheta/dy via Q-tensor
                Qxx=nx.^2-1/2;
                Qxy=nx.*ny;
                [Qxx_x,Qxx_y]=gradient(Qxx);
                [Qxy_x,Qxy_y]=gradient(Qxy);
                dtheta_dx = 0.5.*(Qxx .* Qxy_x - Qxy .* Qxx_x)./(Qxx.^2 + Qxy.^2);
                dtheta_dx(isnan(dtheta_dx)) = 0;
                dtheta_dy = 0.5.*(Qxx .* Qxy_y - Qxy .* Qxx_y)./(Qxx.^2 + Qxy.^2);
                dtheta_dy(isnan(dtheta_dy)) = 0;
                %% calc winding number field
                windingMap = (conv2(dtheta_dy, deltaY, 'same') + conv2(dtheta_dx, deltaX, 'same'))./(2*pi);
                ECMDefects.WindingMap = windingMap;
                %figure, imagesc(windingMap); colorbar
                %% calc defect properties for current loop-integral size
                [positiveCentroids, negativeCentroids, positiveCharges,positiveNematicParameter, negativeCharges, positiveAngle,negativeAngle, negativeNematicParameter] =...
                    DefectAnalysis.calcDefectPositionAngleAndCharge(windingMap, nx, ny, defectCharge, chargeDelta, S,nematicParameterThreshold, nanMap, currentRadiusLineIntegral, ring_kernel);
                %count all defects at current length scale
                if ~isempty(positiveCentroids)
                    numberPosDefectsAtLengthScale(loopInegralSizeIndex) = numberPosDefectsAtLengthScale(loopInegralSizeIndex) + size(positiveCentroids, 1);
                end
                if ~isempty(negativeCentroids)
                    numberNegDefectsAtLengthScale(loopInegralSizeIndex) = numberNegDefectsAtLengthScale(loopInegralSizeIndex) + size(negativeCentroids, 1);
                end
                %%before adding defects from different lengthscale throw
                %%out very nearby defects
                if ~isempty(positiveCentroids_all)
                    throwOutPositive = false(size(positiveCentroids, 1), 1);
                    for k = 1:size(positiveCentroids, 1)
                        if any(mergeDefectsLengthScale^2 > (positiveCentroids(k, 1) - positiveCentroids_all(:, 1)).^2 + (positiveCentroids(k, 2) - positiveCentroids_all(:, 2)).^2)
                            throwOutPositive(k) = true;
                        end
                    end
                    positiveCentroids_all           = [positiveCentroids_all; positiveCentroids(~throwOutPositive, :)];
                    positiveCharges_all             = [positiveCharges_all; positiveCharges(~throwOutPositive)];
                    positiveAngle_all               = [positiveAngle_all; positiveAngle(~throwOutPositive)];
                    positiveNematicParameter_all    = [positiveNematicParameter_all; positiveNematicParameter(~throwOutPositive)];
                else
                    positiveCentroids_all           = [positiveCentroids_all; positiveCentroids];
                    positiveCharges_all             = [positiveCharges_all; positiveCharges];
                    positiveAngle_all               = [positiveAngle_all; positiveAngle];
                    positiveNematicParameter_all    = [positiveNematicParameter_all; positiveNematicParameter];
                    
                end
                
                
                if ~isempty(negativeCentroids_all)
                    throwOutNegative = false(size(negativeCentroids, 1), 1);
                    for k = 1:size(negativeCentroids, 1)
                        if any(mergeDefectsLengthScale^2 > (negativeCentroids(k, 1) - negativeCentroids_all(:, 1)).^2 + (negativeCentroids(k, 2) - negativeCentroids_all(:, 2)).^2)
                            throwOutNegative(k) = true;
                        end
                    end
                    negativeCentroids_all           = [negativeCentroids_all; negativeCentroids(~throwOutNegative, :)];
                    negativeCharges_all             = [negativeCharges_all; negativeCharges(~throwOutNegative)];
                    negativeAngle_all               = [negativeAngle_all; negativeAngle(~throwOutNegative)];
                    negativeNematicParameter_all    = [negativeNematicParameter_all; negativeNematicParameter(~throwOutNegative)];
                else
                    negativeCharges_all             = [negativeCharges_all; negativeCharges];
                    negativeCentroids_all           = [negativeCentroids_all; negativeCentroids];
                    negativeAngle_all               = [negativeAngle_all; negativeAngle];
                    negativeNematicParameter_all    = [negativeNematicParameter_all; negativeNematicParameter];
                end
                
            end
            ECMDefects.positiveCentroids_all                = positiveCentroids_all;
            ECMDefects.negativeCentroids_all                = negativeCentroids_all;
            ECMDefects.positiveCharges_all                  = positiveCharges_all;
            ECMDefects.negativeCharges_all                  = negativeCharges_all;
            ECMDefects.positiveAngle_all                    = positiveAngle_all;
            ECMDefects.negativeAngle_all                    = negativeAngle_all;
            ECMDefects.numberPosDefectsAtLengthScale        = numberPosDefectsAtLengthScale;
            ECMDefects.numberNegDefectsAtLengthScale        = numberNegDefectsAtLengthScale;
            ECMDefects.positiveNematicParameter_all         = positiveNematicParameter_all;
            ECMDefects.negativeNematicParameter_all         = negativeNematicParameter_all;
            
        end
        
        
        
        %------------------------------------------------------------------
        function nanMap = makeNanMap(thetaMap, orientations, nanRegionSizeToAllow)
            %% make NaN-Map in order to delete defects in vicinity
            nanMap = isnan(thetaMap); %use to delete defects in the affected vicinity
            linearInd = sub2ind(size(thetaMap), orientations(:, 2), orientations(:, 1));
            helperMap = zeros(size(nanMap)); helperMap(linearInd) = 1; helperMap = ~helperMap;
            helperMap = bwareaopen(helperMap, nanRegionSizeToAllow);
            nanMap = bwareaopen(nanMap, nanRegionSizeToAllow);
            nanMap = imdilate(nanMap | helperMap, strel('disk', 10));
            
        end
        
        
        
        %------------------------------------------------------------------
        function [deltaX, deltaY, ring_kernel] = makeRingKernel(radiusLineIntegral)
            disky_large = strel('disk',radiusLineIntegral+1);
            disky_large = disky_large.Neighborhood;
            disky_small = strel('disk',radiusLineIntegral);
            disky_small = disky_small.Neighborhood;
            disky_small = padarray(disky_small, (size(disky_large) - size(disky_small))/2);
            ring_kernel= disky_large - disky_small;
            
            edgeLength = length(ring_kernel);
            [X,Y]=meshgrid((1:edgeLength)-edgeLength/2,(1:edgeLength)-edgeLength/2);
            r =(X.^2+Y.^2).^(0.5);
            X(~ring_kernel)=0;
            Y(~ring_kernel)=0;
            deltaX = Y./r;
            deltaY = -X./r;
            
        end
        
         function [dx, dy, square_kernel] = makeSquareKernel(radiusLineIntegral)
             square_kernel = zeros(radiusLineIntegral+1, radiusLineIntegral+1);
             square_kernel(2:(end-1), 2) = 1;
             square_kernel(2:(end-1), end-1) = 1;
             square_kernel(2, 2:(end-1)) = 1;
             square_kernel(end-1, 2:(end-1)) = 1;
             
             [XX, YY] = ind2sub(size(square_kernel), find(square_kernel));
             dx = zeros(size(square_kernel));
             dy = zeros(size(square_kernel));

             dy(2:(end-1), 2)       = -1;
             dy(2:(end-1), end-1)   =  1;
             dx(2, 2:(end-1))       = -1;
             dx(end-1, 2:(end-1))   =  1;
             dy = -dy;
         end
        
        
        %------------------------------------------------------------------
        function [positiveCentroids, negativeCentroids, positiveCharges,positiveNematicParameter, negativeCharges, positiveAngle,negativeAngle, negativeNematicParameter] =...
                calcDefectPositionAngleAndCharge(windingMap, nx, ny, defectCharge, chargeDelta, S,nematicParameterThreshold, nanMap, radiusLineIntegral, ring_kernel)
            %% get defect positions
            mapNegative = windingMap < (-defectCharge+chargeDelta) & windingMap > (-defectCharge - chargeDelta) & S < nematicParameterThreshold & ~nanMap;
            mapPositive = windingMap < (defectCharge+chargeDelta) & windingMap > (defectCharge - chargeDelta) & S < nematicParameterThreshold & ~nanMap;
            %figure, imshow(mapNegative);  figure, imshow(mapPositive);
            %merge very close-by defects (4 pixels --> length scale of kernel)
            structureElement = strel('disk', 2);
            
            mapNegative = imdilate(mapNegative, structureElement);
            mapPositive = imdilate(mapPositive, structureElement);
            mapNegative = imclose(mapNegative, structureElement);
            mapPositive = imclose(mapPositive, structureElement);
            
            CC_neg = bwconncomp(mapNegative);
            CC_pos = bwconncomp(mapPositive);
            
            positiveProps = regionprops(CC_pos);
            positiveCentroids = cat(1, positiveProps.Centroid);
            
            negativeProps = regionprops(CC_neg);
            negativeCentroids = cat(1, negativeProps.Centroid);
            %% calc defect charges
            negativeCentroids = floor(negativeCentroids);
            positiveCentroids = floor(positiveCentroids);
            if ~isempty(negativeCentroids)
                linearIndex_negative        = sub2ind(size(mapNegative), negativeCentroids(:, 2), negativeCentroids(:, 1));
                negativeCharges             = windingMap(linearIndex_negative);
                negativeNematicParameter    = S(linearIndex_negative);
            else
                negativeCharges             = [];
                negativeNematicParameter    = [];
            end
            if ~isempty(positiveCentroids)
                linearIndex_positive            = sub2ind(size(mapNegative), positiveCentroids(:, 2), positiveCentroids(:, 1));
                positiveCharges                 = windingMap(linearIndex_positive);
                positiveNematicParameter        = S(linearIndex_positive);
            else
                positiveCharges             = [];
                positiveNematicParameter    = [];
            end
            %% find defect orientations --> general idea from norton et al.
            positiveAngle = DefectAnalysis.findDefectOrientation(positiveCentroids,ring_kernel,radiusLineIntegral,nx,ny, +1);
            negativeAngle = DefectAnalysis.findDefectOrientation(negativeCentroids,ring_kernel,radiusLineIntegral,nx,ny, -1);
            
            %% throw out defects with charges not in specified range
            posThrowOut                                 = positiveCharges > (defectCharge + chargeDelta) | positiveCharges < (defectCharge - chargeDelta);
            positiveCentroids(posThrowOut, :)           = [];
            positiveAngle(posThrowOut)                  = [];
            positiveCharges(posThrowOut, :)             = [];
            positiveNematicParameter(posThrowOut, :)    = [];
            negThrowOut                                 = negativeCharges > (-defectCharge + chargeDelta) | negativeCharges < (-defectCharge - chargeDelta);
            negativeCentroids(negThrowOut, :)           = [];
            negativeAngle(negThrowOut)                  = [];
            negativeCharges(negThrowOut, :)             = [];
            negativeNematicParameter(negThrowOut, :)    = [];
            
            
        end
        
        
        
        
        %------------------------------------------------------------------
        function phi = findDefectOrientation(centroids,ring,radiusLineIntegral,nx,ny, defectSign)
            %% calc gradients of Q-tensor
            Qxx=nx.^2-1/2;
            Qxy=nx.*ny;
            
            [dxQxx,dyQxx]=gradient(Qxx);
            [dxQxy,dyQxy]=gradient(Qxy);
            
            
            % get coordinates of line integral
            [x1,y1]=find(ring==1);
            x1=x1-(ceil(radiusLineIntegral)+2);
            y1=y1-(ceil(radiusLineIntegral)+2);
            
            numberOfDefects=size(centroids,1);
            phi=zeros(numberOfDefects,1);
            
            [H,W]=size(nx);
            
            for i=1:numberOfDefects
                
                defectX=centroids(i,1);
                defectY=centroids(i,2);
                
                ringAroundDefectX=defectX+x1;
                ringAroundDefectY=defectY+y1;
                
                % handle defects near boundary
                ringAroundDefectX(ringAroundDefectX>W)=W;
                ringAroundDefectX(ringAroundDefectX<1)=1;
                ringAroundDefectY(ringAroundDefectY>H)=H;
                ringAroundDefectY(ringAroundDefectY<1)=1;
                
                
                angleMethod = 3;
                switch angleMethod
                    case 3
                        %% THIRD Method: Vromans, Giomi "Orientational properties of nematic disclinations" --> good results
                        linearIndex_ringAroundDefect = sub2ind(size(nx), ringAroundDefectY, ringAroundDefectX);
                        
                        if defectSign == 1
                            k = 1/2;
                            numerator       = nanmean(dxQxy(linearIndex_ringAroundDefect) - dyQxx(linearIndex_ringAroundDefect));
                            denominator     = nanmean(dxQxx(linearIndex_ringAroundDefect) + dyQxy(linearIndex_ringAroundDefect));
                            phi(i) = k/(1-k)*atan(numerator/denominator);
                        elseif defectSign == -1
                            k = -1/2;
                            numerator       = nanmean(-dxQxy(linearIndex_ringAroundDefect) - dyQxx(linearIndex_ringAroundDefect));
                            denominator     = nanmean(dxQxx(linearIndex_ringAroundDefect) - dyQxy(linearIndex_ringAroundDefect));
                            current_phi = k/(1-k)*atan(numerator/denominator);
                            %take only angles in interval [0, 2/3*pi] for negative
                            %-1/2 defects --> later
                            %                             anglesInDegree = rad2deg([current_phi, current_phi + 2/3*pi, current_phi + 4/3*pi]);
                            %                             anglesInDegree = mod(anglesInDegree, 360);
                            %                             [minVals, ~] = min(anglesInDegree);
                            phi(i) = current_phi;
                        else
                            disp('Defect sign error')
                        end
                        
                end
            end
        end
        
        function plotDefectsOnNativeImage(orgImage, orientations, positiveCentroids, negativeCentroids, positiveAngle, negativeAngle)
            figure, imshow(orgImage)
            vecs = orientations;
            steps = 5;
            hold on,
            directorFieldPlot = quiver(vecs(1:steps:end, 1), vecs(1:steps:end, 2), vecs(1:steps:end, 3),-vecs(1:steps:end, 4), 0.7, "Color", "black");
            directorFieldPlot.ShowArrowHead = false;
            
            defect_arrow_scale = 5;
            
            if ~isempty(negativeCentroids)
                hold on, plot(positiveCentroids(:, 1),positiveCentroids(:, 2), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', [1 0 0])
                centroids_p = positiveCentroids;
                quiver(centroids_p(:,1),centroids_p(:,2),cos(positiveAngle),sin(positiveAngle),defect_arrow_scale,'LineWidth',4, 'Color', [0.2 0.9 0])
                
            end
            if ~isempty(positiveCentroids)
                hold on, plot(negativeCentroids(:, 1),negativeCentroids(:, 2),  'bo', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 1])
                centroids_m = negativeCentroids;
                %plot defect arrows
                negativeAngle = negativeAngle + pi/3;
                quiver(centroids_m(:,1),centroids_m(:,2),cos(negativeAngle),sin(negativeAngle),defect_arrow_scale,'.c','LineWidth',4)
                quiver(centroids_m(:,1),centroids_m(:,2),cos(negativeAngle+2*pi/3),sin(negativeAngle+2*pi/3),defect_arrow_scale,'.c','LineWidth',4)
                quiver(centroids_m(:,1),centroids_m(:,2),cos(negativeAngle+4*pi/3),sin(negativeAngle+4*pi/3),defect_arrow_scale,'.c','LineWidth',4)
            end
            
        end
        
    end
    
end
