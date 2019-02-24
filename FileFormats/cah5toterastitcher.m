classdef  cah5toterastitcher < handle
    %%
    % This class is used to convert a hdf5 dataset to a tiff file
    %
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/04/22
    
    %%
    properties(SetAccess = private, GetAccess = public)
        pathstr;
        h5file;
        flipaxis;
        htiles;
        wtiles;
        dtiles;
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
        background;
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cah5toterastitcher(flipaxis)
            %% Constructor
            % This method requires no argument
            
            %% 
            assert(isa(flipaxis, 'logical'));
            [fullmpath,~,~] = fileparts(mfilename('fullpath'));
            addpath(genpath(fullmpath))
            
            [filename, path] = uigetfile('*.h5', ...
                'Pick the unsticthed data file');
            if filename ~= 0
                obj.pathstr = path;
                obj.flipaxis = flipaxis;
                obj.h5file = cah5file(path, filename, 'r');
                
                obj.htiles = numberofobjects(obj.h5file, '/', 'g');
                if obj.htiles > 0
                    obj.wtiles = numberofobjects(obj.h5file, '/R00', 'g');
                    if obj.wtiles > 0
                        nlast = numberofobjects(obj.h5file, ...
                            sprintf('/R%02d',obj.htiles - 1), 'g');
                        if nlast ~= obj.wtiles
                            obj.htiles = obj.htiles - 1;
                        end
                        nsubstacks = numberofobjects(obj.h5file, ...
                            '/R00/C00', 'd');
                        if nsubstacks > 0
                            dsetname = '/R00/C00/D0000';
                            dimsfirst = ...
                                fliplr(obj.h5file.getdatasetsize(dsetname));
                            dsetname = ...
                                sprintf('/R00/C00/D%04d',64*(nsubstacks-1));
                            dimslast = ...
                                fliplr(obj.h5file.getdatasetsize(dsetname));
                            obj.dtiles = dimsfirst(3) * (nsubstacks - 1) ...
                                + dimslast(3);
                            obj.background = imread([path '/background.tif']);
                        else
                            error(['The selected file is empty or is not' ...
                                'structured as required /Rxx/Cyy/Dzz']);
                        end
                    else
                        error(['The selected file is empty or is not' ...
                            'structured as required /Rxx/Cyy/Dzz']);
                    end
                else
                    error(['The selected file is empty or is not' ...
                        'structured as required /Rxx/Cyy/Dzz']);
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
        
        %%
        function delete(obj)
            %% Destructor
            % This method is called when the object is being destroyed
            % We want to make sure the file is closed
            
            delete(obj.h5file);	% Close the hdf5 file.
        end
        
        %%
        function convert(obj)
            
            progress = 0;
            h = waitbar(0, '0','Name','Creating tiff files ...',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            mapsize = [obj.htiles obj.wtiles];
            ntiles = obj.htiles * obj.wtiles;
            
            tileindex = 1;
            
            dsetname = '/R00/C00/D0000';
            dims = fliplr(obj.h5file.getdatasetsize(dsetname));
            
            while(tileindex <= ntiles )
                % Check for Cancel button press
                if getappdata(h,'canceling')
                    break
                end
                waitbar(progress, h, sprintf('%.2f',progress));
                progress = progress+1/ntiles;
                
                [col, row] = ind2sub(fliplr(mapsize), tileindex);
                
                if obj.flipaxis
                    %y = (obj.wtiles - col) * dims(1);
                    %x = (obj.htiles - row) * dims(2);
                    y = (obj.wtiles - col);
                    x = (obj.htiles - row);
                else
                    %y = (obj.wtiles - col) * dims(1);
                    %x = (obj.htiles - row) * dims(2);
                    y = (obj.wtiles - col);
                    x = (obj.htiles - row);
                end
                
                subfolderName = sprintf('tiff/%06d/%06d_%06d', y, y, x);
                
                subfolder = [obj.pathstr subfolderName];
                
                if ~exist(subfolder, 'dir')
                    mkdir(subfolder);
                end
                
                filename = [subfolder ...
                    sprintf('/%06d_%06d_000000.tif', y, x)];
                
                depth = 0;
                while depth < obj.dtiles
                    if obj.flipaxis
                        datasetName = sprintf('/R%02d/C%02d/D%04d', ...
                            obj.wtiles - col,obj.htiles - row, depth);
                    else
                        datasetName = sprintf('/R%02d/C%02d/D%04d', ...
                            row-1,col-1, depth);
                    end
                    
                    currenttile = readdataset(obj.h5file, datasetName);
                    currenttile = currenttile - repmat(obj.background, [1 1 size(currenttile,3)]);
                    currenttile = max(currenttile, [], 3);
                    
                    if depth == 0
                        if isa(currenttile, 'uint16')
                            bitspersample = 16;
                        elseif isa(currenttile, 'uint8')
                            bitspersample = 8;
                        else
                            error('Only uint8 and uint16 data are supported');
                        end
                        outfile = cabigtiffile(filename, size(currenttile,1),...
                            size(currenttile,2), -1, -1, bitspersample);
                    end
                    
                    outfile.appendstack(currenttile);
                    depth = depth + 64;
                end
                delete(outfile);
                tileindex = tileindex + 1;
                
            end
            delete(h);
        end 
    end
    
end