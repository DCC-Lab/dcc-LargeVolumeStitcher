classdef cabigtiffile < ImageAdapter
    %cabigtiffile - A basic image adapter to read/write Big TIFF files.
    %
    
    %%
    properties(GetAccess = public, SetAccess = private)        
        fullfilename;        
        tiffobj;
        tilelength;
        tilewidth;
        bitspersample;
        empty;
    end
    
    %%
    methods
        
        function obj = cabigtiffile(fname, imagelength, imagewidth, ...
                tilelength, tilewidth, bitspersample)
            %% Constructor
            % tilelength set to -1 for a non tiled file
            % tilewidth set to -1 for a non tiled file
            %%
            validateattributes(fname,       {'char'},   {'row'});
            validateattributes(imagelength, {'numeric'},{'scalar'});
            validateattributes(imagewidth,  {'numeric'},{'scalar'});
            validateattributes(tilelength,  {'numeric'},{'scalar'});
            validateattributes(tilewidth,   {'numeric'},{'scalar'});
            
            if(mod(tilelength,16)~=0 || mod(tilewidth,16)~=0) ...
                    && tilelength ~= -1 && tilewidth ~= -1
                error('cabigtiffile:invalidtilesize',...
                    'Tile size must be a multiple of 16');
            end
            
            obj.fullfilename   = fname;
            obj.ImageSize  = [imagelength, imagewidth, 1];
            obj.tilelength = tilelength;
            obj.tilewidth  = tilewidth;
            obj.bitspersample = bitspersample;
            obj.empty = true;
            
            % Create the Tiff object.
            obj.tiffobj = Tiff(obj.fullfilename, 'w8');
            obj.setTags();          
        end
        
        function delete(obj)
            %% Destructor
            % This method is called when the object is being destroyed
            % We want to make sure the file is closed
            
            %%
            obj.close();
        end
        
        function setTags(obj)
            % Setup the tiff file properties
            obj.tiffobj.setTag('ImageLength',   obj.ImageSize(1));
            obj.tiffobj.setTag('ImageWidth',    obj.ImageSize(2));
            if obj.tilelength ~= -1 && obj.tilewidth ~= -1
                obj.tiffobj.setTag('TileLength',    obj.tilelength);
                obj.tiffobj.setTag('TileWidth',     obj.tilewidth);
            end
            obj.tiffobj.setTag('Photometric',   ...
                Tiff.Photometric.MinIsBlack);
            obj.tiffobj.setTag('BitsPerSample', obj.bitspersample);
            obj.tiffobj.setTag('SampleFormat',  Tiff.SampleFormat.UInt);
            obj.tiffobj.setTag('SamplesPerPixel', 1);
            obj.tiffobj.setTag('PlanarConfiguration', ...
                Tiff.PlanarConfiguration.Chunky); 
        end
        
        
        function [] = write2Darray(obj, offset, wdata)
            %% Write a array to the file.
            
            %%
            assert(mod(offset(1), obj.tilelength) == 0);
            assert(mod(offset(2), obj.tilewidth) == 0);
            if mod(offset(1), obj.tilelength)~=0 ...
                    || mod(offset(2), obj.tilewidth)~=0
                error('cabigtiffile:invalidoffset',...
                    'Offset must match a tile start');
            end
            
            address = offset;
            tilesize = [obj.tilelength obj.tilewidth];
            dataptr = [1:tilesize(1); 1:tilesize(2)];
            
            
            while true
                % offset to a tile number.
                tilenumber = obj.tiffobj.computeTile(address(1:2));
                
                % Loop until all data are written
                folderIdx = 0;
                while folderIdx < size(wdata,3)
                    % If wdata is less than tile size, this function it will
                    % pads with 0s.
                    dirNum = offset + folderIdx;
                    obj.tiffobj.setDirectory(dirNum)
                    data = wdata(dataptr(1,:), dataptr(2,:),:);
                    obj.tiffobj.writeEncodedTile(tilenumber, data);
                    folderIdx = folderIdx + 1;
                end
                
                %
                dataptr(1,:) = dataptr(1,:) + tilesize(1);
                address(1) = address(1) + tilesize(1);
                if dataptr(1,end) > size(wdata,1)
                    dataptr(1,:) = 1:tilesize(1);
                    dataptr(2,:) = dataptr(2,:) + tilesize(2);
                    address(2) = address(2);
                    if dataptr(2,end) > size(wdata,2)
                        break;
                    end
                end
                
            end
            
        end
        
        %%
        function [] = appendstack(obj, wdata)
            %% Append a substack of images to the file
            % 
            %
            
            %%
            assert(size(wdata,1) == obj.ImageSize(1) && ...
                size(wdata,2) == obj.ImageSize(2));
            sliceindex = 1;
            slicecount = size(wdata,3);
            
            if obj.empty
                obj.tiffobj.write(wdata(:,:,1));
                sliceindex = sliceindex + 1;
                obj.empty = false;
            end
            
            while sliceindex <= slicecount
                obj.tiffobj.writeDirectory();
                obj.setTags();
                obj.tiffobj.write(wdata(:,:,sliceindex));
                sliceindex = sliceindex + 1;
            end
        end
        
        
        %%
        function rdata = readslice(obj, sliceIndex)
            %% Read a slice
            % Call this method to read a dataset.The method takes one
            % argument:
            % sliceIndex - 
            
            %%
            obj.tiffobj.setDirectory(sliceIndex);
            rdata = obj.tiffobj.read();
        end
        
        %%
        function addslice(obj, wdata)
            %% create a new directory and add a slice
            % Call this method to create a new directory and write data to
            % it. The method takes two arguments:
            % wdata - array to write.
            
            %%
            assert(size(wdata,1) == obj.ImageSize(1) && ...
                size(wdata,2) == obj.ImageSize(2));
            if ~obj.empty
                obj.tiffobj.writeDirectory();
                obj.setTags();
            else
                obj.empty = false;
            end
            obj.tiffobj.write(wdata);
            
        end
        
        
        function close(obj)
            % Close the tiff file
            obj.tiffobj.close();
        end
        
        function readRegion(~)
        end
        
    end    
end