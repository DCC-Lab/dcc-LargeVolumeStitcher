classdef  cah5tobigdataviewer < handle
    %%
    % Class to convert our hdf5 file to ImageJ Big Data Wiewer
    % format
    %
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/06/27
    
    %%
    properties(SetAccess = private, GetAccess = public)
        h5file;
        name;
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cah5tobigdataviewer(varargin)
            %% Create a new cah5tobigdataviewer object            
            %%
            if nargin > 0
                assert(ischar(varargin{1}));
                [path, filename,~] = fileparts(varargin{1});
                filename = [filename '.h5'];
            else
                [filename, path] = uigetfile('*.h5','Pick the a hdf5 file');
            end
            if filename ~= 0
                obj.h5file = cah5file(path, filename, 'a');
                obj.name = [path, filename];
            else
                warndlg('No file selected');
            end
        end
        
        %%
        function delete(obj)
            %% Called when objected is being deleted
            % This method is called when the object is being destroyed
            % We want to make sure the file is closed
            
            delete(obj.h5file);	% Close the hdf5 file.
        end
        
        %%
        function convert(obj)
            %% Perform the conversion to imagej big data viewer format
            % Need a data set following the path /t00000/s00/0/cells
            
            %%
            dsetname = '/t00000/s00/0/cells';
            dims = fliplr(obj.h5file.getdatasetsize(dsetname));
            
            try
                resolutions = obj.h5file.readdataset('/s00/resolutions');
            catch
                resolutions = [1 1 1];
            end
            
            try
                subdivisions = obj.h5file.readdataset('/s00/subdivisions');
            catch
                subdivisions = [32 32 16];
            end
            h = waitbar(0,'0', 'Name', 'Writing the pyramidal file ...');
            
            for sliceindex = 1:16:dims(3)
                offset = [0 0 sliceindex-1];
                remainingSlices = dims(3) - sliceindex + 1;
                if remainingSlices < 16
                    block = [dims(1:2) remainingSlices];
                else
                    block = [dims(1:2) 16];
                end
                slicedata = obj.h5file.readblock(dsetname, offset, block);
                res = 1;
                while min(size(slicedata, 1), size(slicedata,2)) > 512
                    %slicedata = slicedata(1:2:end,1:2:end);
                    slicedata = imresize(slicedata, 0.5, 'bilinear');
                    offset = offset/2;
                    %if mod(sliceindex-1,2^res) == 0
                        newDstName = sprintf('t00000/s00/%d/cells', res);
                        appendtodataset(obj.h5file, newDstName, offset, slicedata);
                    %end
                    res = res + 1;
                end
                waitbar(sliceindex/dims(3),h, sprintf('%.2f',sliceindex/dims(3)));
                clear slicedata;
            end
            resolutions = single(resolutions *  2.^(0:res-1));
            write2ddataset(obj.h5file, '/s00/resolutions', resolutions);
            
            subdivisions = single(subdivisions * ones(1, res));
            write2ddataset(obj.h5file, '/s00/subdivisions', subdivisions);
            
            close(h);
        end
        
        function N = binning(~, M, p)
            [m,n]=size(M); %M is the original matrix
            
            M=sum( reshape(M,p,[]) ,1 );
            M=reshape(M,m/p,[]).'; %Note transpose
            
            M=sum( reshape(M,p,[]) ,1);
            N=reshape(M,n/p,[])'; %Note transpose
        end
   
        
    end
    
end