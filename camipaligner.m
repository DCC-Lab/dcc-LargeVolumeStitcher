classdef camipaligner < handle
    %%
    % This class is used to compute pairwise translation
    % between a image/stack with its south and east neighbour.
    % Stacks may not exceed 128 images
    % The translations are obtained by:
    % 1. For each stack, compute maximum projection intensity(MIP) along
    %    all directions.
    % 2. For each direction, compute translation using a combination of phase
    %    correlation and normalized cross-correlation
    % 3. For each axis, select the translation associated with the smallest
    %    euclidian square distance
    
    %%
    properties(SetAccess = public, GetAccess = public)
        reference;  % Reference image/stack
        template;   % template image/stack
    end
    
    %%
    properties(SetAccess = private, GetAccess = public)
        % translation of the template image/stack relative to reference
        translation;
        % correlation coefficients associatated to translation
        correlationcoef;
        % Quadratic error
        sqerror;
    end
    
    %%
    methods(Access = public)
        %%
        function obj = camipaligner(varargin)
            %% Create a new camipaligner object
            % Arguments
            % arg 1. Reference stack/image
            % arg 2. Template stack/image
            % arg 3. Ratio of the overlap area
            %
            % obj = camipaligner(ref, template, overlap) create
            %       the object and compute translations
            
            %%
            assert(nargin > 3);
            assert(ischar(varargin{4}));
            
            obj.reference = varargin{1};
            obj.template = varargin{2};
            
            [obj.translation, obj.correlationcoef, obj.sqerror] = ...
                obj.mipalign(varargin{4}, varargin{3});
        end
        
        %%
        function free(obj)
            %% Release the memory used by the object
            % This method is used to free the memory used by the
            % object. All the properties will be set to empty
            
            %%
            obj.reference = [];
            obj.template     = [];
            
            obj.translation = [];
            obj.correlationcoef = [];
        end
        
    end
    
    %%
    methods(Access = private)     
        %%
        function [u,v,c,err] = pcncc(obj, subject, template)
            %% Calculate translation between 2 images using PC and NCC
            %
            % Combine phase correlation nad fast normalized
            % cross-correlation to find translation between 2 images of
            % same size.
            %
            % returns
            % * u : rows shift
            % * v : cols shift
            % * c : correlation associated to this shift
            % * er: rms error after shift adjustment
            % see Reference below for more details
            %
            % Y. Yu and H. Peng, ?Automated high speed stitching of large
            % 3d microscopic images,?
            % in IEEE 2011 International Symposium on Biomedical Imaging:
            % From Nano to Macro (ISBI?2011),
            % Los Alamitos,CA, USA,2011,pp.238?241,IEEE Computer Society.
            %
            
            %%
            dimensions = [ndims(subject) ndims(template)];
            assert(...
                (dimensions(1) == dimensions(2))...
                && ((dimensions(1) == 2) || (dimensions(1) == 3 && size(subject,3)))...
                );
            % Initialize features for I
            if size(subject,3) == 3
                subject = rgb2gray(subject);
                template = rgb2gray(template);
            end
            % Image size
            rows = size(subject,1);
            cols = size(subject,2);
            % Apodization
            %win = window(@chebwin, size(subject,1)) * window(@chebwin, size(subject,2))';
            win = hanning(size(subject,1)) * hanning(size(subject,2))';
            subject = double(subject).* win;
            template = double(template).* win;
            % Fourier tranform
            sfft = fft2(double(subject));
            tfft = fft2(double(template));
            % Phase correlation
            cps = sfft.*conj(tfft);
            ncps = cps./abs(cps);
            pc = abs(ifft2(ncps));
            pc(isnan(pc)) = 0;
            % find peaks and keep the 4 highest
            npeaks = 8;
            pcLocs = find(imregionalmax(pc));
            [~, pindex] = sort(pc(pcLocs), 'descend');
            npindex = numel(pindex);
            if npindex >= npeaks
                pcLocs = pcLocs(pindex(1:npeaks));
            else
                pcLocs = pcLocs(pindex);
                npeaks = npindex;
            end
            
            % Normalized cross correlation
            nPixels = numel(subject);
            sSum = sum(subject(:));
            tSum = sum(template(:));
            s2Sum = sum(subject(:).^2);
            t2Sum = sum(template(:).^2);
            
            numer = abs(ifft2(cps)) - tSum*sSum/nPixels;
            denom= sqrt((t2Sum - tSum^2/nPixels)*(s2Sum - sSum^2/nPixels));
            ncc = numer/denom;
            
            % Find exact value of the peak using NCC in a 5x5 region
            % surrounding peaks found by PC
            nccLocs = zeros(npeaks,2);
            nccPeaks = zeros(npeaks,1);
            for n = 1:npeaks
                pcLoc = [rem(pcLocs(n), rows) floor(pcLocs(n)/rows) + 1];
                if(pcLoc(1) == 0)
                    pcLoc(1) = rows;
                    pcLoc(2) = pcLoc(2) - 1;
                end
                
                row = [pcLoc(1)-9, pcLoc(1)+9];
                col = [pcLoc(2)-9, pcLoc(2)+9];
                if row(1) < 1, row(1) = 1; end
                if row(2) > size(ncc,1), row(2) = size(ncc,1);end
                if col(1) < 1, col(1) = 1; end
                if col(2) > size(ncc,2), col(2) = size(ncc,2);end
                subNCC = ncc(row(1):row(2), col(1):col(2));
                nccPeaks(n) = max(subNCC(:));
                nccLocs(n,:) = pcLoc;
            end;
           
            % return the highest
            [c, ccIdx] = max(nccPeaks);
            u = nccLocs(ccIdx,1);
            v = nccLocs(ccIdx,2);
            
            if u > fix(rows/2), u = u - rows - 1;
            else u = u - 1;
            end
            
            if v > fix(cols/2), v = v - cols - 1;
            else v = v - 1;
            end
            
            err = obj.rmserror(subject, template, u, v);
            
        end
        
        %%
        function [u,v,c,err] = xcc(obj, subject, template)
            %% Calculate translation between 2 images using cross correlatio
            %
            % use cross-correlation to find translation between 2 images of
            % same size.
            %
            % returns
            % * u : rows shift
            % * v : cols shift
            % * c : correlation associated to this shift
            % * er: rms error after shift adjustment
            
            %%
            dimensions = [ndims(subject) ndims(template)];
            assert(...
                (dimensions(1) == dimensions(2))...
                && ((dimensions(1) == 2) || (dimensions(1) == 3 && size(subject,3)))...
                );
            % Initialize features for I
            if size(subject,3) == 3
                subject = rgb2gray(subject);
                template = rgb2gray(template);
            end
            % Image size
            rows = size(subject,1);
            cols = size(subject,2);
            
            %apodization
            win = hanning(size(subject,1)) * hanning(size(subject,2))';
            subject = single(subject).* win;
            template = single(template).* win;
            
            % Fourier transform
            subject =subject - mean2(subject);
            template = template - mean2(template);
            sfft = fft2(subject);
            tfft = fft2(template);
            cc = ifft2(sfft.*conj(tfft))/(norm(subject) * norm(template)); 
            
            [max1,loc1] = max(cc);
            [~,loc2] = max(max1);
            u=loc1(loc2);
            v=loc2;
            c=cc(u,v);
            
            if u > fix(rows/2), u = u - rows - 1;
            else u = u - 1;
            end
            
            if v > fix(cols/2), v = v - cols - 1;
            else v = v - 1;
            end
            
            err = obj.rmserror(subject, template, u, v);
            
        end
        
        %%
        function err = rmserror(~, subject, template, u, v)
            lenu = size(subject, 1);
            lenv = size(subject, 2);
            if u >=  0
                su = u + 1:lenu;
                tu = 1:lenu - u;
            else
                su = 1:lenu + u;
                tu = 1 - u:lenu;
            end
            
            if v >=  0
                sv = [v + 1 lenv];
                tv = [1 lenv - v];
            else
                sv = 1:lenv + v;
                tv = 1 - v:lenv;
            end
            
            den = 1/rms(rms(subject + template));
            d = subject(su, sv) - template(tu, tv);
            err = rms(d(:))* den;
        end
        
        %%
        function [u, c, err] = mipalign(varargin)
            %% Calculate shift between 2 stacks using MIP and (PC+NCC) or XCC
            %
            % Use maximum intensity projections and cross-correlation to 
            % obtan the shift between 2
            % stacks of images.
            %
            % The stacks must be previously loaded in the object and
            % additionnaly, the method needs:
            % * direction : the relative position of the subject stack
            %   relative to the reference one: 'E' for east, 'S' for south
            % * overlap   : ration of the overlaping region
            %
            % Optionaly, user can choose PCNCC method in place of XCC, by
            % adding 'p' or 'P' as 4th argument.
            % The strategy of PCNCC is detailed in the reference paper below. 
            % Unlike the paper, we are using PNCC (instead of PC) to find
            % the correlation peaks
            %
            % returns
            % * u : shift (width)
            % * v : shift (height)
            % * w : shift (depth)
            % * c : matrix containing correlation associated with the shift
            %       in each direction
            %
            % NOTES: TRIED PARALLEL COMPUTING BUT IT WAS 2X LONGER THAN
            %        SEQUENTIAL
            
            %%
            assert(nargin > 2);
            obj = varargin{1};
            direction = varargin{2};
            overlap = varargin{3};
            if nargin > 3, meth = varargin{4}; else meth = 'x'; end
            assert(ischar(meth));
            assert(~isempty(obj.reference));
            assert(upper(direction) == 'E' || upper(direction) == 'S');
            rows = ceil(size(obj.reference,1) * overlap);
            cols = ceil(size(obj.reference,2) * overlap);
            %frames = size(obj.reference,3);
            
            % Select the overlaping region
            assert(~isempty(obj.template));
            switch upper(direction)
                case 'E'
                    ref = obj.reference(:,end-cols+1:end, :);
                    scene = obj.template(:,1:cols, :);
                case 'S'
                    ref = obj.reference(end-rows+1:end,:, :);
                    scene = obj.template(1:rows,:, :);
                otherwise
            end
            
            % Calculate maximum intensities projections for ref. stack
            %             mip23Ref = reshape(max(ref, [], 1), size(ref,2), frames);
            %             mip13Ref = reshape(max(ref, [], 2), size(ref,1), frames);
            mip12Ref = max(ref, [], 3);
            % Calculate maximum intensities projections for scene stack
            %mip23Scene = reshape(max(scene, [], 1), size(ref,2),frames);
            %mip13Scene = reshape(max(scene, [], 2), size(ref,1), frames);
            mip12Scene = max(scene, [], 3);
            
            if upper(meth) == 'P'
                % phase correlation - normalized cross-correlation
                %[col1, depth1, cc23, err23] = pcncc(obj, mip23Ref, mip23Scene);
                %[row1, depth2, cc13, err13] = pcncc(obj, mip13Ref, mip13Scene);
                [row2, col2, cc12, err12] = pcncc(obj, mip12Ref, mip12Scene);
            else
                % Cross correlation
                %[col1, depth1, cc23, err23] = xcc(obj, mip23Ref, mip23Scene);
                %[row1, depth2, cc13, err13] = xcc(obj, mip13Ref, mip13Scene);
                [row2, col2, cc12, err12] = xcc(obj, mip12Ref, mip12Scene);
            end
            
            % Select the displacement corresponding to the highest peak
            %             if frames > 63
            %                 if err13 < err12, u(1) = row1; c(1) = cc13; err(1) = err13;
            %                 else u(1) = row2; c(1) = cc12; err(1) = err12;
            %                 end
            %                 if err23 < err12, u(2) = col1; c(2) = cc23; err(2) = err23;
            %                 else u(2) = col2; c(2) = cc12; err(2) = err12;
            %                 end
            %                 if err23 < err13, u(3) = depth1; c(3) = cc23; err(3) = err23;
            %                 else u(3) = depth2; c(3) = cc13; err(3) = err13;
            %                 end
            %             else
            %                 u = [row2 col2 0];
            %                 c = [cc12 cc12 0];
            %                 err = [err12 err12 0];
            %             end
            u = [row2 col2 0];
            c = [cc12 cc12 0];
            err = [err12 err12 0];
            
            switch upper(direction)
                case 'E'
                    u = u - [0 cols 0]; % the real shift includes overlap
                case 'S'
                    u = u - [rows 0 0]; % the real shift includes overlap
                otherwise
            end
            
        end
    end
end