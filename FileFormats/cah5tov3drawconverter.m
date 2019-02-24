classdef  cah5tov3drawconverter < handle
    %%
    % This class is used to convert a hdf5 dataset to a v3d's raw file
    %
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/03/24
    
    %%
    properties(SetAccess = private, GetAccess = public)
        h5file;
        v3drawfile;
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cah5tov3drawconverter()
            %% Constructor
            % This method requires no argument
            
            %% 
            [filename, path] = uigetfile('*.h5','Pick the a hdf5 file');
            if filename ~= 0
                obj.h5file = cah5file(path, filename, 'r');
                [~,name,~] = fileparts(filename);
                outfilename = [name '.v3draw'];
                obj.v3drawfile = cav3drawfile(path, outfilename, 'w');
                
            else
                warndlg('No file selected');
            end
        end
        
        %%
        function delete(obj)
            %% Destructor
            % This method is called when the object is being destroyed
            % We want to make sure the file is closed
            
            delete(obj.h5file);	% Close the hdf5 file.
            delete(obj.v3drawfile);	% Close the v3d file.
        end
        
        %%
        function convert(obj, dsetname)
            dims = fliplr(obj.h5file.getdatasetsize(dsetname));
            info.bytesPerPixel = obj.h5file.getbytespersample(dsetname);
            
            info.width = dims(2);
            info.height = dims(1);
            info.depth = dims(3);
            info.timeSetups = 1;
            
            h = waitbar(0,'CONVERTING TO V3DRAW');
            
            obj.v3drawfile.writeinfo(info);
            for sliceindex = 1:dims(3)
                offset = [0 0 sliceindex-1];
                block = [dims(1:2) 1];
                slicedata = obj.h5file.readblock(dsetname, offset, ...
                    block);
                obj.v3drawfile.writeslice(sliceindex, slicedata);
                waitbar(sliceindex/dims(3),h);
            end
            close(h);
        end
   
        
    end
    
end