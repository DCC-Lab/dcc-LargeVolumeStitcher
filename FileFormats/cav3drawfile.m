classdef  cav3drawfile < handle
    %%
    % This class is used to open/create a v3d's raw file
    %
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/03/24
    
    %%
    properties(SetAccess = private, GetAccess = public)
        fullfilename; % The full file specification.
        info;
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
        fileID; % the file ID of the successfully opened/created file
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cav3drawfile(path, filename, option)
            %% Constructor
            % This method requires 3 arguments
            % path - path to the location of the file
            % filename - name of the file with extension
            % option - to indicate whether the file is to be created ('w')
            %          or opened ('r') for read
            
            %%
            assert(nargin == 3);
            
            obj.fullfilename = fullfile(path, filename);
            assert(ischar(obj.fullfilename)==1);
            
            assert(option == 'r' || option == 'w');
            
            if option == 'r'
                obj.fileID = fopen(obj.fullfilename, 'r');
                obj.info = obj.getinfo();
            elseif option == 'w'
                obj.fileID = fopen(obj.fullfilename, 'w');
            end
            assert(obj.fileID ~= -1);
        end
        
        function delete(obj)
            %% Destructor
            % This method is called when the object is being destroyed
            % We want to make sure the file is closed
            
            fclose(obj.fileID);	% Close the file.
        end
        
        %%
        function info = getinfo(obj)
            %% Read V3D's raw file header
            
            % Go to the beginning of the file.
            fseek(obj.fileID, 0, 'bof');
            
            %
            info.format = fread(obj.fileID,24, '*char')';
            % Read endianness
            info.endianness = fread(obj.fileID,1, '*char');
            
            % Read bytes per pixel
            info.bytesPerPixel = fread(obj.fileID,1, '*int16');
            
            % Read in the dimensions for the different directions.
            info.width = fread(obj.fileID, 1, '*int32');
            info.height = fread(obj.fileID, 1, '*int32');
            info.depth = fread(obj.fileID, 1, '*int32');
            info.timeSetups = fread(obj.fileID, 1, '*int32');
            nbrPixels = info.width * info.height;
            info.pixelsPerSlice = int32(nbrPixels);
            info.bytesPerSlice = int32(info.bytesPerPixel) * nbrPixels;
        end
        
        %%
        function writeinfo(obj, info)
            %% Write V3D's raw file header
            % info must be a structure with the following fields
            % bytesPerPixel - 
            % width - 
            % height -
            % depth -
            % timeSetups - 
            
            assert(isstruct(info));
            assert(isfield(info, 'bytesPerPixel'));
            assert(isfield(info, 'width'));
            assert(isfield(info, 'height'));
            assert(isfield(info, 'depth'));
            assert(isfield(info, 'timeSetups'));
            
            info.format = 'raw_image_stack_by_hpeng';
            info.endianness = 'L';
            
            % Go to the beginning of the file.
            fseek(obj.fileID, 0, 'bof');
            % Write format
            fwrite(obj.fileID, info.format', '*char');
            % write endianless
            fwrite(obj.fileID, info.endianness, '*char');            
            % Write bytes per pixel
            fwrite(obj.fileID, info.bytesPerPixel, '*int16');             
            % Write in the dimensions for the different directions.
            dims = [info.width; info.height; info.depth; info.timeSetups];
            fwrite(obj.fileID, dims, '*int32');
            
            info.bytesPerSlice = info.bytesPerPixel * info.width * info.height;
            obj.info = info;
        end
        
        %%
        function rdata = readslice(obj, sliceIndex)
            %% Read a slice
            % Call this method to read a dataset.The method takes one
            % argument:
            % sliceIndex - 
            
            %%
            offset = 43 + obj.info.bytesPerSlice * (sliceIndex - 1);
            fseek(obj.fileID, offset ,'bof');
            switch obj.info.bytesPerPixel
                case 1
                    rdata = fread(obj.fileID, obj.info.pixelsPerSlice,...
                        '*uint8');
                case 2
                    rdata = fread(obj.fileID, obj.info.pixelsPerSlice,...
                        '*uint16'); 
                case 4
                    rdata = fread(obj.fileID, obj.info.pixelsPerSlice,...
                        '*uint32');
                otherwise
            end
        end
        
        %%
        function writeslice(obj, sliceIndex, wdata)
            %% Write a slice
            % Call this method to create a new dataset and write data to
            % it. The method takes two arguments:
            % sliceIndex - full path form the root group to the dataset
            % wdata - array to write.
            
            %%
            assert(size(wdata,1) == obj.info.height && ...
                size(wdata,2) == obj.info.width);
            
            offset = 43 + obj.info.bytesPerSlice * (sliceIndex - 1);
            fseek(obj.fileID, offset ,'bof');
            switch class(wdata)
                case 'uint8'
                    fwrite(obj.fileID, uint16(wdata'), '*uint8'); 
                case 'uint16'
                    fwrite(obj.fileID, uint16(wdata'), '*uint16'); 
                case 'uint32'
                    fwrite(obj.fileID, uint16(wdata'), '*uint32'); 
                case 'double'
                    fwrite(obj.fileID, wdata', '*double'); 
                otherwise
            end
            
        end
        
    end
    
end