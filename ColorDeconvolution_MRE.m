classdef ColorDeconvolution_MRE
    %following the Paper of Macenko et al.
    methods(Static)
        function OD = calcOD(Value, Maximum)
            %ODs scale linearly with amount of stain
            OD = max(0, -log10(Value/Maximum));
        end
        function ODtable =  calcODtable(MaximumBackground)
            %create OD lookup table for performace
            format long
            ODtable = zeros(255, 1);
            Range = 0.0001:1/255:1;
            for i = 1:255
                ODtable(i) = ColorDeconvolution_MRE.calcOD(Range(i), MaximumBackground); %zeros at the end of ODtable correspond to MaximumBackground
            end
        end
        
        function Covariance = calcCovariance(x,y)
            x_mean = mean(x);
            y_mean = mean(y);
            Covariance = 0;
            for i = 1:length(x)
                x_err = x(i) - x_mean;
                y_err = y(i) - y_mean;
                Covariance = Covariance + x_err*y_err/(length(x));
            end
        end
        
        function [stain1, stain2, Residuum] = estimateStainsInImage(im, MaximumBackground, ignorePercentage, HE, grayThresh)
            
            ignorePercentage = ignorePercentage/100;
            im = double(im);
            im = reshape(im, size(im, 1)*size(im,2), 3); %shape into rgb columns
            im(im == 0) = 1; % im=1 --> looked in ODtable corresponds to OD of 0.0001 Intensity
            ODtable = ColorDeconvolution_MRE.calcODtable(MaximumBackground);
            OD_RGB_Columns = ODtable(im); %converst RGB columns to ODs
            format long
            redOD = OD_RGB_Columns(:,1);
            greenOD = OD_RGB_Columns(:,2);
            blueOD = OD_RGB_Columns(:,3);
            if HE == true
                %% HE specific filters
%               HE_Filter_Eosin     = ~(im(:,1) > im(:,2) & im(:,3) > im(:,2));
               HE_Filter_Hemat     = im(:,3) > im(:,2) & im(:,3) > im(:,1);
%               HE_Filter_Bloody    = im(:,1) > im(:,3) & im(:,3) > im(:,2) & im(:,1) > im(:,2);
                HE_Filter = (im(:,1) > im(:,2) | im(:,3) > im(:,2)) | HE_Filter_Hemat;
                %% general filters for colors
                MaximumBackground_sq = MaximumBackground*MaximumBackground;
                Pxl_sq               = (im(:,1).*im(:,1)./(255^2) + im(:,2).*im(:,2)./(255^2) + im(:,3).*im(:,3)./(255^2));
                sq3                  = 1/sqrt(3);
                GrayFilter = 3*MaximumBackground_sq < Pxl_sq | Pxl_sq <= 0 | (im(:,1)*sq3 + im(:,2)*sq3 + im(:,3)*sq3)./(Pxl_sq).^(0.5) <= 3*grayThresh;
                 %% too white or too black
                 maxIm = max(im(:));
                WhiteOrBlack = (im(:,1) < 10 & im(:,2) < 10 & im(:,3) < 10) | (im(:,1) > 0.98*maxIm & im(:,2) > 0.98*maxIm & im(:,3) > 0.98*maxIm);
                SortOut =  GrayFilter | HE_Filter | WhiteOrBlack; %| HE_Filter_Bloody;
                
                redOD(SortOut) = [];
                greenOD(SortOut) = [];
                blueOD(SortOut) = [];
                OD_RGB_Columns(SortOut, :) = [];
            else %do normal filtering
                %% general filters for colors
                maxIm = max(im(:));
                MaximumBackground_sq = MaximumBackground*MaximumBackground;
                Pxl_sq               = (im(:,1).*im(:,1)./(255^2) + im(:,2).*im(:,2)./(255^2) + im(:,3).*im(:,3)./(255^2));
                sq3                  = 1/sqrt(3);
                GrayFilter = MaximumBackground_sq < Pxl_sq | Pxl_sq == 0 | (im(:,1)*sq3 + im(:,2)*sq3 + im(:,3)*sq3)./(Pxl_sq).^(0.5) <= grayThresh |  im(:,2) > .5.*(im(:,1) + im(:,3));
                WhiteOrBlack = (im(:,1) < 10 & im(:,2) < 10 & im(:,3) < 10) | (im(:,1) > 0.95*maxIm & im(:,2) > 0.95*maxIm & im(:,3) > 0.95*maxIm);
                SortOut = GrayFilter | WhiteOrBlack;
                
                redOD(SortOut) = [];
                greenOD(SortOut) = [];
                blueOD(SortOut) = [];
                OD_RGB_Columns(SortOut, :) = [];
            end
            CovarianceMatrix = zeros(3,3);
            CovarianceMatrix(1,1) = ColorDeconvolution_MRE.calcCovariance(redOD, redOD);
            CovarianceMatrix(2,2) = ColorDeconvolution_MRE.calcCovariance(greenOD, greenOD);
            CovarianceMatrix(3,3) = ColorDeconvolution_MRE.calcCovariance(blueOD, blueOD);
            CovarianceMatrix(1,2) = ColorDeconvolution_MRE.calcCovariance(redOD, greenOD);
            CovarianceMatrix(1,3) = ColorDeconvolution_MRE.calcCovariance(redOD, blueOD);
            CovarianceMatrix(2,3) = ColorDeconvolution_MRE.calcCovariance(greenOD, blueOD);
            CovarianceMatrix(3,1) = CovarianceMatrix(1,3);
            CovarianceMatrix(3,2) = CovarianceMatrix(2,3);
            CovarianceMatrix(2,1) = CovarianceMatrix(1,2);
            [V,D] = eig(CovarianceMatrix); %eigenwert problem
            [~, ind] = sort(unique(D(D~=0))); %two biggest eigenvalues
            eigenV1 = V(:,ind(3));
            eigenV2 = V(:,ind(2));
            [theta,~] = cart2pol(redOD*eigenV1(1) + greenOD*eigenV1(2) + blueOD*eigenV1(3)...
                ,redOD*eigenV2(1) + greenOD*eigenV2(2) + blueOD*eigenV2(3));
            [~, Indices] = sort(theta);
            %find two angles which span the color space
            LowerSpan = Indices(max(1, ceil(ignorePercentage*length(theta))));
            UpperSpan = Indices(min(length(theta)-1, floor((1-ignorePercentage)*length(theta))));
            %extract stain vectors
            stain1 = OD_RGB_Columns(LowerSpan, :);
            stain1 = stain1./norm(stain1);
            stain2 = OD_RGB_Columns(UpperSpan, :);
            stain2 = stain2./norm(stain2);
            Residuum = cross(stain1, stain2);
            Residuum = Residuum./norm(Residuum);
        end
        
        function [Channel_1, Channel_2, Residuum_Channel, stain1, stain2] = calcColorTransform(im, MaximumBackground, IgnorePercentiles, HE, grayThresh, image_Downsample, preCalculatedStains)
            %% INPUTS Main function
            %           im:                uint8 image
            %           MaxBackground:     UpperThresh
            %           IgnorePercentiles: Ignore Percentiles of angles spanning the color space
            %           HE:                true if HE image, false if not
            %           grayThresh:        thresh for gray pixels where
            %           red=green=blue  < graythresh
            %           image_Downsample:      between 0 and 1, is the
            %           downscaling factor for speeding up the estimation
            %           of the stain vectors (1 is original size)
            % preCalculatedStains: if false stains are estimated else give
            % 3*3 matrix with colums
            
            %% get image-specific stain vectors
            downSample_im = imresize(im , image_Downsample);
            if preCalculatedStains  == false
                [stain1, stain2, Residuum] = ColorDeconvolution_MRE.estimateStainsInImage(downSample_im,MaximumBackground, IgnorePercentiles, HE, grayThresh);
            else
               stain1 = preCalculatedStains( 1, :);
               stain2 = preCalculatedStains(2, :);
               Residuum = preCalculatedStains(3, :);
            end
            im = double(im);
            im (im == 0) = 1;
            sz1 = size(im,1);
            sz2 = size(im,2);
            im = reshape(im, sz1*sz2, 3);
            ODtable = ColorDeconvolution_MRE.calcODtable(MaximumBackground);
            OD_RGB_Columns = ODtable(im); %converst RGB columns to ODs
            StainConcentration = OD_RGB_Columns*inv([stain1; stain2; Residuum]); %estimate how much stain is in each pixel
            StainConcentration = reshape(StainConcentration, sz1, sz2, 3);
            if HE == true
                [~, max_ind] = max([stain1(1) stain2(1)]);
                Channel_1 =  StainConcentration(:,:,max_ind); %--> Hematoxylin
                if max_ind == 1
                    Channel_2 =  StainConcentration(:,:,2);%--> Eosin
                elseif max_ind == 2
                    Channel_2 =  StainConcentration(:,:,1);%--> Eosin
                end
                Residuum_Channel = StainConcentration(:,:,3);%-->Residuum
                
            else %one has to find out later which stain is which
                Channel_1 = StainConcentration(:,:,1);
                Channel_2 = StainConcentration(:,:,2);
                Residuum_Channel = StainConcentration(:,:,3);
            end
        end
    end
end