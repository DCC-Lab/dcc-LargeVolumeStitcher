classdef cavolumestitcher < handle
    %%
    % This class is used to compute pairwise translation between 3d stacks
    % of images and stitch them together.
    % The file contanining these images must be a hdf5 file that respect
    % the following dataset nomenclature and organization:
    % /Rxx/Cyy/Dzzzz where xx is the row coordinate of the tile, yy is the
    % column coordinate and zzzz the depth coordinate. 
    % Large stacks are splitted into substacks of 64 images each. Dzzzz
    % indicate the position of the substack in the large stack
    %
    %% Usage
    % 1.Create an object
    %   obj = cavolumestitcher(flipaxis)
    %   will prompt to choose the file contanining unstitched images.
    %   flipaxis is true or false to indicate that axis must be flipped,
    %   Origin is supposed to be located in the top left conner. If it's
    %   not the case for your data, use the flipaxis property.
    %   At this stage, the output file will be created based on the name of
    %   the 'input file name' + 'stitched' + '.h5'
    %
    % 2. Perform the stitching operation
    %    obj.stitch();
    %    This can take many seconds/minutes to complete depending on the
    %    size of the input file. Note that the time consuming operation
    %    here is writing output to disk, so using SSD drive will help (will
    %    be 4x faster).
    %    It is possible to cancel and restart the stitching without
    %    deleting the object. Use the button on the progress bar or 
    %   
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/03/21
    
    %%
    properties(SetAccess = private, GetAccess = public)
        % HDF5 file containing non-stitched tiles
        inputfile;
        % HDF5 file where stitched slices will be written
        outputfile;
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
        % Wait bar handle
        progressbar;
        % tiles absolute positions
        abspositions;
        % open tiles
        inmemorytiles;
        % absolute translations
        globalshifts;
        % 3d matrix containing row, column & depth translation of each tile
        pairwiseshifts;
        % 3d matrix containing correlation associated to the translations
        ccpairwise;
        % 
        pairwisesqerr;
        % 3d matrix
        tilesadjancymat;
        % number of columns of tiles
        wtiles;
        % number of rows of tiles
        htiles;
        % number of substacks per tile
        dtiles;
        % true if axis are to be flipped, false otherwise
        flipaxis;
        %
        imrows;
        %
        imcols;
        %
        background;
        % Histogram
        %histogram;
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cavolumestitcher(flipaxis)
            %% Constructor
            %
            %
            
            %%
            assert(isa(flipaxis, 'logical'));
            [fullmpath,~,~] = fileparts(mfilename('fullpath'));
            addpath(genpath(fullmpath))
            
            [filename, path] = uigetfile('*.h5', ...
                'Pick the unsticthed data file');
            if filename ~= 0
                obj.flipaxis = flipaxis;
                obj.inputfile = cah5file(path, filename, 'r');
                
                obj.htiles = numberofobjects(obj.inputfile, '/', 'g');
                if obj.htiles > 0
                    obj.wtiles = numberofobjects(obj.inputfile, '/R00', 'g');
                    if obj.wtiles > 0
                        nlast = numberofobjects(obj.inputfile, ...
                            sprintf('/R%02d',obj.htiles - 1), 'g');
                        if nlast ~= obj.wtiles
                            obj.htiles = obj.htiles - 1;
                        end
                        obj.dtiles = numberofobjects(obj.inputfile, ...
                            '/R00/C00', 'd');
                        if obj.dtiles > 0
                            [~,name,~] = fileparts(filename);
                            outfilename = [name 'stiched.h5'];
                            outpath = [path '/processed'];
                            if 7~=exist(outpath,'dir')
                                mkdir(path,'processed');
                            end
                            obj.outputfile = cah5file(outpath, outfilename, 'w');
                        else
                            error(['The selected file is empty or is not' ...
                                'structured as required /Rxx/Cyy/Dzz']);
                        end
                    else
                        error(['The selected file is empty or is not' ...
                            'structured as required /Rxx/Cyy/Dzzzz']);
                    end
                else
                    error(['The selected file is empty or is not' ...
                        'structured as required /Rxx/Cyy/Dzzzz']);
                end
                
                if obj.flipaxis
                    tmp = obj.htiles;
                    obj.htiles = obj.wtiles;
                    obj.wtiles = tmp;
                end
                
            else
                error('Not enough arguments');
            end
        end
        
        function delete(obj)
            %% Called when the object is being destroyed
            
            %%
            if ishandle(obj.progressbar)
                delete(obj.progressbar);
            end
        end
        
        function stitch(obj)
            %% Perform the stitching operation
            tic
            success = obj.computepairwisetranslations();
           if success == 1
                obj.optimizemap();
                obj.computeabstranslation(1);
                obj.computeabstranslation(2);
                %obj.computeabstranslation(3);
                obj.tileabsoluteposition();
                obj.buildmap();
                %close files
                outfilename = obj.outputfile.fullfilename;
                obj.inputfile.delete();
                obj.outputfile.delete();
                % pyramids
                pyramidFile = cah5tobigdataviewer(outfilename);
                pyramidFile.convert();
                pyramidFile.delete();
            end
            toc
        end
        
    end
    
    %%
    methods(Access = private)
        %%
%         function success = computegainandoffset(obj)
%             % 1. read first tile
%             % 2. neighbours
%             % 3. check if in memory
%             % 4. read missing neighbour
%             % 5. mippcncc
%             % 6. clear current
%             % 7. load next
%             % 8. go to 2.
%             
%             success = true;
%             progress = 0;
%             if ishandle(obj.progressbar) 
%                 delete(obj.progressbar);
%             end
%             obj.progressbar = waitbar(0.0, '0','Name',...
%                 'Computing gain and offset ...',...
%                 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%             
%             mapsize = [obj.htiles obj.wtiles];
%             ntiles = obj.htiles * obj.wtiles;
%              
%             tileindex = 1;
%             depth = 64;
%             %kk = 1;
%             cumsumarray = [];
%             nelt = 0;
%                       
%             while ~isempty(obj.inputfile)
%                 % Check for Cancel button press
%                 if getappdata(obj.progressbar,'canceling')
%                     success = false;
%                     break
%                 end
%                 waitbar(progress, obj.progressbar, sprintf('D%d - %.5f', depth, ...
%                     progress));
%                 progress = progress+1/ntiles;
%                                 
%                 [row, col] = ind2sub(mapsize, tileindex);
%                 if obj.flipaxis
%                     datasetName = sprintf('/R%02d/C%02d/D%04d', ...
%                         obj.wtiles - col,obj.htiles - row, depth);
%                 else
%                     datasetName = sprintf('/R%02d/C%02d/D%04d', ...
%                         row-1,col-1, depth);
%                 end
%                 
%                 currenttile = readdataset(obj.inputfile, datasetName);
%                 if ~isempty(cumsumarray)
%                     cumsumarray = cumsumarray + sum(single(currenttile),3);
%                 else
%                     cumsumarray = sum(single(currenttile),3);
%                 end
%                 nelt = nelt + size(currenttile,3);
%                 
%                 tileindex = tileindex + 1;
%                 
%                 if tileindex > ntiles
%                     depth = depth + 64;
%                     kk = kk + 1;
%                     tileindex = 1;
%                     progress = 0;
%                     if kk > obj.dtiles
%                         break;
%                     end
%                 end
%             end
%             
%             if isa(currenttile, 'uint8')
%                 obj.background = uint8(medfilt2(cumsumarray/nelt, [5,5]));
%             elseif isa(currenttile, 'uint16')
%                 obj.background = uint16(medfilt2(cumsumarray/nelt, [5,5]));
%             end;
%             
%             delete(obj.progressbar);
%         end
        %%
        function success = computepairwisetranslations(obj)
            % 1. read first tile
            % 2. neighbours
            % 3. check if in memory
            % 4. read missing neighbour
            % 5. mippcncc
            % 6. clear current
            % 7. load next
            % 8. go to 2.
            
            success = true;
            progress = 0;
            if ishandle(obj.progressbar) 
                delete(obj.progressbar);
            end
            obj.progressbar = waitbar(0.0, '0','Name',...
                'Computing pairwise translations ...',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            mapsize = [obj.htiles obj.wtiles];
            ntiles = obj.htiles * obj.wtiles;
            
            obj.pairwiseshifts  = zeros(ntiles, ntiles, 3);
            obj.pairwisesqerr = ones(ntiles, ntiles, 3);
            obj.tilesadjancymat = zeros(ntiles, ntiles, 3);
            %obj.histogram = [];
            
            tileindex = 1;
            depth = 64;
            %kk = 1;
            cumsumarray = [];
            nelt = 0;
            
            try
                tilesOverlap = readattribute(obj.inputfile, [], 'tilesOverlap');
            catch
                prompt = {'Enter tiles overlap [0.0 - 1.0]:'};
                dlg_title = 'Tiles overlap missing from the input file';
                num_lines = 1;
                defaultans = {'0.1','hsv'};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                tilesOverlap = str2double(answer);
            end
            
            while ~isempty(obj.inputfile)
                % Check for Cancel button press
                if getappdata(obj.progressbar,'canceling')
                    success = false;
                    break
                end
                waitbar(progress, obj.progressbar, sprintf('D%d - %.5f', depth, ...
                    progress));
                progress = progress+1/ntiles;
                                
                [row, col] = ind2sub(mapsize, tileindex);
                if obj.flipaxis
                    datasetName = sprintf('/R%02d/C%02d/D%04d', ...
                        obj.wtiles - col,obj.htiles - row, depth);
                else
                    datasetName = sprintf('/R%02d/C%02d/D%04d', ...
                        row-1,col-1, depth);
                end
                
                currenttile = readdataset(obj.inputfile, datasetName);
                if ~isempty(obj.background)
                    currenttile = currenttile - ...
                        repmat(obj.background, [1 1 size(currenttile,3)]);
                else
                    if ~isempty(cumsumarray)
                        cumsumarray = cumsumarray + sum(single(currenttile),3);
                    else
                        cumsumarray = sum(single(currenttile),3);
                    end
                    nelt = nelt + size(currenttile,3);
                end
                
                %                 if isa(currenttile, 'uint16') || isa(currenttile, 'int16')
                %                     nbins = 65536;
                %                 else
                %                     nbins = 256;
                %                 end
                %
                %                if isempty(obj.histogram)
                %                    [obj.histogram, ~]= histcounts(currenttile,0:nbins-1);
                %                else
                %                    [theHist, ~] = histcounts(currenttile,0:nbins-1);
                %                    obj.histogram = obj.histogram + theHist;
                %                end
                
                obj.imrows = size(currenttile,1);
                obj.imcols = size(currenttile,2);
                
                jj = sub2ind(mapsize, row, col);
                
                if row ~= 1
                    northstitcher = camipaligner(obj.inmemorytiles{row-1}, ...
                        currenttile, tilesOverlap, 'S');
                    ii = sub2ind(mapsize,row-1,col);
                    obj.pairwiseshifts(ii,jj,:)  = northstitcher.translation;
                    obj.pairwisesqerr(ii,jj,:) = northstitcher.sqerror; %northstitcher.correlationcoef;
                end
                
                if col ~= 1
                    weststitcher = camipaligner(obj.inmemorytiles{row}, ...
                        currenttile, tilesOverlap, 'E');
                    ii = sub2ind(mapsize,row,col-1);
                    obj.pairwiseshifts(ii,jj,:)  = weststitcher.translation;
                    obj.pairwisesqerr(ii,jj,:) = weststitcher.sqerror; %weststitcher.correlationcoef;
                end
                obj.inmemorytiles{row} = currenttile;
                tileindex = tileindex + 1;
                
                if tileindex > ntiles
                    %depth = depth + 64;
                    %kk = kk + 1;
                    %tileindex = 1;
                    %progress = 0;
                    %if kk > obj.dtiles
                    break;
                    %end
                end
            end
            
            if isa(currenttile, 'uint8')
                obj.background = uint8(medfilt2(cumsumarray/nelt, [5,5]));
            elseif isa(currenttile, 'uint16')
                obj.background = uint16(medfilt2(cumsumarray/nelt, [5,5]));
            end;
            
            %if success
            %    [obj.pairwisesqerr, minIdx] = min(obj.pairwisesqerr, [], 3);
            %    obj.pairwiseshifts = squeeze(obj.pairwiseshifts(minIdx));
            %    obj.pairwisesqerr = squeeze(obj.pairwisesqerr);
            %end
            delete(obj.progressbar);
        end
        
        %%
        function optimizemap(obj)           
            
            %             A = sparse(obj.pairwisesqerr(:,:,1));
            %             B = A + A';
            %             G = graph(B);
            %             T = minspantree(G);
            %             obj.tilesadjancymat(:,:,1) = adjacency(T);
            %
            %             A = sparse(obj.pairwisesqerr(:,:,2));
            %             B = A + A';
            %             G = graph(B);
            %             T = minspantree(G);
            %             obj.tilesadjancymat(:,:,2) = adjacency(T);
            %
            %             A = sparse(obj.pairwisesqerr(:,:,3));
            %             B = A + A';
            %             G = graph(B);
            %             T = minspantree(G);
            %             obj.tilesadjancymat(:,:,3) = adjacency(T);

            ntiles = obj.wtiles * obj.htiles;
            for nn = 1: ntiles-1
                [~, idx] = min(obj.pairwisesqerr(nn,:,1));
                obj.tilesadjancymat(nn,idx,1) = 1;
                [~, idx] = min(obj.pairwisesqerr(nn,:,2));
                obj.tilesadjancymat(nn,idx,2) = 1;
                [~, idx] = min(obj.pairwisesqerr(nn,:,3));
                obj.tilesadjancymat(nn,idx,3) = 1;
            end
            
        end
        
        %%
        function buildmap(obj)
            
            ntiles = obj.wtiles * obj.htiles;
            mapsize = [obj.htiles obj.wtiles];
            
            tileindex = 1;
            depth = 0;
            
            progress = 0;
            if ishandle(obj.progressbar) 
                delete(obj.progressbar);
            end
            obj.progressbar = waitbar(progress/ntiles, '0.00','Name','Building map...',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            try
                voxelSize = readattribute(obj.inputfile, [], 'voxelSize');
            catch
                voxelSize = single([1 1 1]);
            end
            
            % work on histogram;
            % find maximum in the histogram and remove bins upper than it
            % 256 bins
            % keep bins with counts greater than 250
            % This is the maximum for normalization
            %             cumhist = cumsum(obj.histogram);
            %             cumhist = cumhist/max(cumhist);
            %             histMaxCount = find(cumhist==1,1);
            %             obj.histogram(histMaxCount:end) = [];
            %             oldHistLen = length(obj.histogram);
            %             newHistBinLen = ceil(oldHistLen/256);
            %             newHistLen = 256*newHistBinLen;
            %             obj.histogram(oldHistLen+1:newHistLen) = 0;
            %             obj.histogram = obj.downsamplebin(obj.histogram,2,newHistBinLen);
            %             cumhist = cumsum((obj.histogram > 250));
            %             cumhist = cumhist/max(cumhist);
            %             maxNorm = find(cumhist==1,1) * newHistBinLen;
            
            while ~isempty(obj.inputfile)
                % Check for Cancel button press
                if getappdata(obj.progressbar,'canceling')
                    break
                end
                waitbar(progress, obj.progressbar, sprintf('D%d - %.2f', depth, progress));
                progress = progress + 1/ntiles;
                                
                [row, col] = ind2sub(mapsize, tileindex);
                if obj.flipaxis
                    datasetName = sprintf('/R%02d/C%02d/D%04d', ...
                        obj.wtiles - col,obj.htiles - row, depth);
                else
                    datasetName = sprintf('/R%02d/C%02d/D%04d', ...
                        row-1,col-1, depth);
                end
                
                currenttile = readdataset(obj.inputfile, datasetName);
                if ~isempty(obj.background)
                    bgnd = double(obj.background);
                    gain = repmat(mean2(bgnd)./bgnd, ...
                        [1 1 size(currenttile,3)]);
                    offset = repmat(bgnd, [1 1 size(currenttile,3)]);
                    %currenttile = double(currenttile);
                    %offset = repmat(mean(currenttile,3), [1 1 size(currenttile,3)]);
                    %currenttile = uint16((currenttile - offset).*gain);
                end
                
                %currenttile = uint16(double(currenttile) * 65535/maxNorm);
                
                if col ~= 1
                    west = obj.inmemorytiles{row};
                    [westrowidx, westcolidx, rowindex, colindex] = ...
                        obj.overlappingarea(row, col-1, row, col);
                    if ~isempty(westrowidx) && ~isempty(westcolidx)
                        fused = obj.alphablend(...
                            west(westrowidx,westcolidx,:), ...
                            currenttile(rowindex,colindex,:), 'y');
                        currenttile(rowindex,colindex,:) = fused;
                        west(westrowidx,westcolidx,:) = fused;
                        clear fused;
                    end                  
                end
                
                if row ~= 1
                    north = obj.inmemorytiles{row-1};
                    [northridx, northcidx, rowindex, colindex] = ...
                        obj.overlappingarea(row-1, col, row, col);
                    if ~isempty(northridx) && ~isempty(northcidx)
                        fused = obj.alphablend(...
                            north(northridx,northcidx,:), ...
                            currenttile(rowindex,colindex,:), 'x');
                        currenttile(rowindex,colindex,:) = fused;
                        north(northridx,northcidx,:) = fused;
                        
                        obj.inmemorytiles{row - 1} = north;
                        clear fused;
                    end
                end
                
                obj.inmemorytiles{row} = currenttile;
                
                if exist('west', 'var')  
                    absrow = obj.abspositions(row, col-1, 1);
                    abscol = obj.abspositions(row, col-1, 2);
                    %%
                    dstName = '/t00000/s00/0/cells';
                    appendtodataset(obj.outputfile, dstName, ...
                        [absrow-1, abscol-1, depth], west);
                end
                
                tileindex = tileindex + 1;
                
                %%
                if tileindex > ntiles
                    colidx = obj.wtiles;
                    for rowidx = 1:obj.htiles
                        tile = obj.inmemorytiles{rowidx};
                        absrow = obj.abspositions(rowidx, colidx, 1);
                        abscol = obj.abspositions(rowidx, colidx, 2);
                        %%
                        dstName = '/t00000/s00/0/cells';
                        appendtodataset(obj.outputfile, dstName, ...
                            [absrow-1, abscol-1, depth], tile);
                    end
                    depth = depth + 64;
                    tileindex = 1;
                    clear west;
                    progress = 0;
                    if bitsra(depth,6) >= obj.dtiles
                        break;
                    end
                end
            end
            
            %createattribute(obj.outputfile, '/t00000/s00/0/cells', ...
            %    'element_size_um', voxelSize(2));
            
            write2ddataset(obj.outputfile, '/s00/resolutions', fliplr(voxelSize)');
            subdivisions = single([32 32 16]);
            write2ddataset(obj.outputfile, '/s00/subdivisions', subdivisions');
            delete(obj.progressbar);
        end
        
        %%
        function computeabstranslation(obj, axe)
            assert(axe < 4 && axe > 0);
            % Find out vertices in the spanning tree
            [vertices(:,1), vertices(:,2)] = find(obj.tilesadjancymat(:,:,axe));
            % first tile is the central one
            %
            node = 1;
            myTree{1} = []; % intialize the tree
            treeIdx = 1;    % no branches for the moment
            ramif = [];     % ramification nodes
             
            while (1)
                % all direct vertices from/to current node
                [childrow, childcol]  = find(vertices==node);
                if ~isempty(childrow)
                    % many vertices from/to current nodes: we have a branch
                    % here
                    if numel(childrow) > 1
                        ramif = [ramif node];
                        N = numel(childrow) - 1;
                        treeLen = length(myTree);
                        % Create a new branch and copy out the path from 
                        % node 1 to the ramification node
                        for nn = 1:N
                            myTree{treeLen+nn} = myTree{treeIdx};
                        end
                    end
                    if childcol(1) == 1, nextcol = 2;
                    else nextcol = 1; end
                    node = vertices(childrow(1),nextcol);
                    k = childrow(1);
                    
                    tsign = 1 - (nextcol == 1) * 2; % translation dir
               
                    % copy the vertex and the corresponding translation sign
                    % in the tree [node row col sign]
                    % node is the node to which translation is applied
                    % [row, col] cordinates of translation in the t matrix
                    % sign indicate the dir translation is applied
                    myTree{treeIdx} = ...
                        [myTree{treeIdx}; [node vertices(k,:) tsign]];
                    vertices(k,:) = [];
                elseif ~isempty(ramif)
                    node = ramif(1);
                    ramif(1) = [];
                    treeIdx = treeIdx + 1;
                else
                    break;
                end
            end
            
            tMatrix = obj.pairwiseshifts(:,:, axe);
            
            for ii = 1:length(myTree)
                vect = myTree{ii};
                if ~isempty(vect)
                    tVect = ...
                        diag(tMatrix(vect(:,2), vect(:,3))) .* vect(:,4);
                    tVect = cumsum(tVect);
                    obj.globalshifts(vect(:,1), axe) = tVect;
                end
            end
        end
        
        %%
        function tileabsoluteposition(obj)
            
            tmat(:,:,1) = reshape(obj.globalshifts(:,1), obj.htiles, obj.wtiles);
            tmat(:,:,2) = reshape(obj.globalshifts(:,2), obj.htiles, obj.wtiles);
            
            firstRow = tmat(1,:,1);
            firstCol = tmat(:,1,2);
            
            refrow = -min(firstRow);
            refcol = -min(firstCol);
            
            ntiles = obj.wtiles * obj.htiles;
            mapsize = [obj.htiles obj.wtiles];
            for tile = 1:ntiles
                [row, col] = ind2sub(mapsize, tile);
                absrow = refrow + (row - 1) * obj.imrows + tmat(row, col, 1) + 1;
                abscol = refcol + (col - 1) * obj.imcols + tmat(row, col, 2) + 1;
                obj.abspositions(row, col,:) = [absrow abscol];
            end
        end
        
        %%
        function [rrows, rcols, trows, tcols] = overlappingarea(obj, ...
                rrow, rcol, trow, tcol)
            %% overlaping region
            % @arguments
            
            rpos = squeeze(obj.abspositions(rrow, rcol, :));
            tpos = squeeze(obj.abspositions(trow, tcol, :));
            
            rrows = rpos(1) : rpos(1) + obj.imrows-1;
            rcols = rpos(2) : rpos(2) + obj.imcols-1;
            
            trows = tpos(1) : tpos(1) + obj.imrows-1;
            tcols = tpos(2) : tpos(2) + obj.imcols-1;
            
            rrindex = ismember(rrows,trows);
            rcindex = ismember(rcols,tcols);
            
            trindex = ismember(trows,rrows);
            tcindex = ismember(tcols,rcols);
            
            rrows = find(rrindex);
            rcols = find(rcindex);
            
            trows = find(trindex);
            tcols = find(tcindex);
        end
        
        %%
        function fused = alphablend(~, ref, temp, blendaxis)
            assert(ischar(blendaxis));
            ref = double(ref);
            temp = double(temp);
            [nx,ny,nz] = size(ref);
            blendaxis = lower(blendaxis);
            if blendaxis == 'y'
                X = ones(nx, 1);
                Y = linspace(0,1,ny);
            else
                X = linspace(0,1,nx)';
                Y = ones(1, ny);
            end
            W = X * Y;
            W = repmat(W, [1 1 nz]);
            fused = uint16(W.*temp + (1-W).*ref);
            %fused = max(ref,temp); 
        end
        
        %%
        function [outv] = downsamplebin(~,invec,dimnum,amount)
            sizebefore = size(invec);
            sizemiddle = [sizebefore(1:dimnum-1),amount,sizebefore(dimnum)/amount,sizebefore(dimnum+1:length(sizebefore))];
            sizeafter = sizebefore;
            sizeafter(dimnum) = sizeafter(dimnum)/amount;
            outv = reshape(sum(reshape(invec,sizemiddle),dimnum),sizeafter);
        end
        
    end
    
end