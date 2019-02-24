classdef  cah5totiffconverter < handle
    %%
    % This class is used to convert a hdf5 dataset to a tiff file
    %
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/04/22
    
    %%
    properties(SetAccess = private, GetAccess = public)
        h5file;
        outfilename;
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cah5totiffconverter()
            %% Constructor
            % This method requires no argument
            
            %% 
            [filename, path] = uigetfile('*.h5','Pick the a hdf5 file');
            if filename ~= 0
                obj.h5file = cah5file(path, filename, 'r');
                [~,name,~] = fileparts(filename);
                obj.outfilename = [path '/' name]; 
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
        end
        
        %%
        function convert(obj, dsetname)
                
            dims = fliplr(obj.h5file.getdatasetsize(dsetname));
            [~,name] = fileparts(dsetname);
            outfile = cabigtiffile([obj.outfilename '_' name '.tif'], ...
                dims(1), dims(2), -1, -1, 16);
            
            h = waitbar(0,'0', 'Name', 'Converting to tif ...');
            
            for sliceindex = 1:1:dims(3)
                offset = [0 0 sliceindex-1];
                block = [dims(1:2) 1];
                slicedata = obj.h5file.readblock(dsetname, offset, ...
                    block);
                %slicedata = max(slicedata, [], 3);
                outfile.addslice(slicedata);
                waitbar(sliceindex/dims(3),h, sprintf('%.2f',...
                    sliceindex/dims(3)));
            end
            close(h);
            delete(outfile);
        end
   
        
    end
    
end