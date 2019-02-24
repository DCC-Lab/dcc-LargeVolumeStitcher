classdef  cah5file < handle
    %%
    % This class is used to open/create a hdf5 file
    %
    % Author: Cleo. Akitegetse
    % Copyright 2016 Cleo. Akitegetse
    % Date: 2016/03/21
    
    %%
    properties(SetAccess = private, GetAccess = public)
        fullfilename; % The full file specification.
    end
    
    %%
    properties(SetAccess = private, GetAccess = private)
        fileID; % the file ID of the successfully opened/created file
    end
    
    %%
    methods(Access = public)
        %%
        function obj = cah5file(path, filename, option)
            %% Create a new file object
            % Arguments:
            %   path     - path to the location of the file
            %   filename - name of the file with extension
            %   option   - to indicate whether the file is to be created ('w')
            %              or opened ('r') for read
            
            %%
            assert(option == 'r' || option == 'w' || option =='a');
            
            obj.fullfilename = fullfile(path, filename);
            
            switch option
                case 'r'
                    obj.fileID = H5F.open(obj.fullfilename,'H5F_ACC_RDONLY',...
                        'H5P_DEFAULT');
                case 'w'
                    fcpl = H5P.create('H5P_FILE_CREATE');
                    fapl = H5P.create('H5P_FILE_ACCESS');
                    obj.fileID = H5F.create(obj.fullfilename, ...
                        'H5F_ACC_TRUNC',fcpl,fapl);
                case 'a'
                    obj.fileID = H5F.open(obj.fullfilename,'H5F_ACC_RDWR',...
                        'H5P_DEFAULT');
            end
            assert(obj.fileID ~= -1);
        end
        
        function delete(obj)
            %% Called when the object is being destroyed
            % We want to make sure the file is closed
            
            %%
            H5F.close(obj.fileID);
        end
        
        %%
        function n = numberofobjects(obj, groupName, type)
            %% Count the number of objects of the specified type in the group
            % Arguments:
            %   groupName - the name of the group containing objects
            %   type      - type of objects
            %               'l' - link, 'g' - group,
            %               'd' - dataset, 't' - type
            
            %%
            plist = 'H5P_DEFAULT';
            gid = H5G.open(obj.fileID, groupName);
            ginfo = H5G.get_info(gid);
            idxtype = 'H5_INDEX_NAME';
            order = 'H5_ITER_DEC';
            n = 0;
            for j = 0:ginfo.nlinks-1
                oid = H5O.open_by_idx(obj.fileID,groupName,idxtype,order,j,plist);
                oinfo = H5O.get_info(oid);
                switch(lower(type))
                    case 'l'
                        if oinfo.type == H5ML.get_constant_value('H5G_LINK')
                            n = n + 1;
                        end
                        
                    case 'g'
                        if oinfo.type == H5ML.get_constant_value('H5G_GROUP')
                            n = n + 1;
                        end
                        
                    case 'd'
                        if oinfo.type == H5ML.get_constant_value('H5G_DATASET')
                            n = n + 1;
                        end
                    case 't'
                        if oinfo.type == H5ML.get_constant_value('H5G_TYPE')
                            n = n + 1;
                        end
                end
                H5O.close(oid);
            end
            H5G.close(gid);
        end
        
        %%
        function rdata = readdataset(obj, datasetName)
            %% Read a dataset
            % Call this method to read a dataset.The method takes one
            % Arguments:
            %    datasetName - full path from the root group to the dataset
            
            %%
            assert(isa(datasetName, 'char'));
            assert(obj.fileID ~= -1);
            dset = H5D.open (obj.fileID, datasetName);
            assert(dset ~= -1);
            rdata = H5D.read(dset,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',...
                'H5P_DEFAULT');
            H5D.close(dset);
        end
        
        %%
        function bdata = readblock(obj, datasetName, offset, block)
            %% Read a block in a dataset
            % Call this method to read a dataset.The method takes three
            % Arguments:
            %   datasetName - full path from the root group to the dataset
            %   offset      - location of the block
            %   block       - size of the block
            
            %%
            assert(isa(datasetName, 'char'));
            assert(obj.fileID ~= -1);
            dset = H5D.open (obj.fileID, datasetName);
            assert(dset ~= -1);
            
            % how the data is to be selected from the file
            filespace = H5D.get_space(dset);
            [nDims, ~] = H5S.get_simple_extent_dims(filespace);
            offset = fliplr(offset);
            block = fliplr(block);
            assert(nDims == numel(offset));
            H5S.select_hyperslab(filespace,'H5S_SELECT_SET',offset,[],...
                [],block);
            % how the data is to be arranged in memory
            memspace = H5S.create_simple(nDims,block,[]);
            % read data
            bdata = H5D.read(dset,'H5ML_DEFAULT',memspace,filespace,...
                'H5P_DEFAULT');
            % Close resources
            H5D.close(dset);
            H5S.close(memspace);
            H5S.close (filespace);
        end
        
        %%
        function writedataset(obj, datasetName, wdata)
            %% Write a 3D array to a dataset
            % Call this method to create a new dataset and write data to
            % it.
            % Arguments:
            %   datasetName - full path from the root group to the dataset
            %   wdata       - array to write to the dataset.
            
            %%
            assert(isa(datasetName, 'char'));
            assert(obj.fileID ~= -1);
            assert(ndims(wdata)<4);
            % Create the data space with unlimited dimensions.
            if ndims(wdata) == 3
                dims = fliplr(size(wdata));
            else
                dims = fliplr([size(wdata) 1]);
            end
            unlim = H5ML.get_constant_value('H5S_UNLIMITED');
            maxdims = [unlim unlim unlim];
            dataspace = H5S.create_simple(3,dims,maxdims);
            % Data type
            switch class(wdata)
                case 'uint16'
                    datatype = H5T.copy('H5T_NATIVE_UINT16');
                case 'uint8'
                    datatype = H5T.copy('H5T_NATIVE_UINT8');
                case 'int16'
                    datatype = H5T.copy('H5T_NATIVE_INT16');
                case 'int8'
                    datatype = H5T.copy('H5T_NATIVE_INT8');
                case 'double'
                    datatype = H5T.copy('H5T_NATIVE_DOUBLE');
                case 'single'
                    datatype = H5T.copy('H5T_NATIVE_FLOAT');
                otherwise
                    datatype = H5T.copy('H5T_NATIVE_INT');
            end
            % Modify dataset creation properties, i.e. enable chunking
            dcplid = H5P.create('H5P_DATASET_CREATE');
            lcplid = H5P.create('H5P_LINK_CREATE');
            H5P.set_create_intermediate_group(lcplid,1);
            daplid = H5P.create('H5P_DATASET_ACCESS');
            H5P.set_chunk(dcplid,fliplr([32 32 16]));
            % Create a new dataset within the file using chunk creation
            % properties.
            dataset = H5D.create(obj.fileID, datasetName, datatype, ...
                dataspace, lcplid, dcplid, daplid);
            assert(dataset ~= -1);
            % Write data to dataset
            H5D.write(dataset,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',...
                'H5P_DEFAULT',wdata);
            % Close resources
            H5D.close(dataset);
            H5P.close(dcplid);
            H5P.close(lcplid);
            H5P.close(daplid);
            H5T.close(datatype);
            H5S.close (dataspace);
        end
        
        %%
        function write2ddataset(obj, datasetName, wdata)
            %% Write a 2D array to a dataset
            % Call this method to create a new dataset and write data to
            % it.
            % Arguments:
            %   datasetName - full path from the root group to the dataset
            %   wdata       - array to write to the dataset.
            
            %%
            assert(isa(datasetName, 'char'));
            assert(obj.fileID ~= -1);
            assert(ismatrix(wdata));
            
            dims = fliplr(size(wdata));
            % Create the data space with unlimited dimensions.
            unlim = H5ML.get_constant_value('H5S_UNLIMITED');
            maxdims = [unlim unlim];
            dataspace = H5S.create_simple(2,dims,maxdims);
            
            % Data type
            switch class(wdata)
                case 'uint16'
                    datatype = H5T.copy('H5T_NATIVE_UINT16');
                case 'uint8'
                    datatype = H5T.copy('H5T_NATIVE_UINT8');
                case 'int16'
                    datatype = H5T.copy('H5T_NATIVE_INT16');
                case 'int8'
                    datatype = H5T.copy('H5T_NATIVE_INT8');
                case 'double'
                    datatype = H5T.copy('H5T_NATIVE_DOUBLE');
                case 'single'
                    datatype = H5T.copy('H5T_NATIVE_FLOAT');
                otherwise
                    datatype = H5T.copy('H5T_NATIVE_INT');
            end
            % Modify dataset creation properties, i.e. enable chunking
            dcplid = H5P.create('H5P_DATASET_CREATE');
            lcplid = H5P.create('H5P_LINK_CREATE');
            H5P.set_create_intermediate_group(lcplid,1);
            daplid = H5P.create('H5P_DATASET_ACCESS');
            H5P.set_chunk(dcplid,fliplr([32 32]));
            
            % Check if the dataset exists already
            dataset = -1;
            try
                dataset = H5D.open(obj.fileID,datasetName);
            catch
            end
            
            
            if dataset ~= -1
                H5D.set_extent(dataset,dims);
            else
                % Create a new dataset within the file using chunk creation
                % properties.
                dataset = H5D.create(obj.fileID, datasetName, datatype, ...
                    dataspace, lcplid, dcplid, daplid);
                assert(dataset ~= -1);
            end
            % Write data to dataset
            H5D.write(dataset,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',...
                'H5P_DEFAULT',wdata);
            % Close resources
            H5D.close(dataset);
            H5P.close(dcplid);
            H5P.close(lcplid);
            H5P.close(daplid);
            H5T.close(datatype);
            H5S.close (dataspace);
        end
        
        %%
        function appendtodataset(obj, datasetName, offset, wdata)
            %% Append a 3D array to an existing dataset
            % Call this method to append data to an existing dataset.
            % Arguments:
            %   datasetName - full path from the root group to the dataset
            %   offset      - offset location of the array in the existing dataset
            %   wdata       - array to append.
            
            %%
            assert(isa(datasetName, 'char'));
            assert(obj.fileID ~= -1);
            assert(ndims(wdata)<4);
            % Open dataset
            dataset = -1;
            try
                dataset = H5D.open(obj.fileID,datasetName);
            catch
            end
            
            if(dataset ~= -1)
                typeId = H5D.get_type(dataset);
                % Data type
                switch class(wdata)
                    case 'uint16'
                        assert(H5T.equal(typeId, 'H5T_NATIVE_UINT16') == 1);
                    case 'uint8'
                        assert(H5T.equal(typeId, 'H5T_NATIVE_UINT8') == 1);
                    case 'int16'
                        assert(H5T.equal(typeId, 'H5T_NATIVE_INT16') == 1);
                    case 'int8'
                        assert(H5T.equal(typeId, 'H5T_NATIVE_INT8') == 1);
                    case 'double'
                        assert(H5T.equal(typeId, 'H5T_NATIVE_DOUBLE') == 1);
                    case 'single'
                        assert(H5T.equal(typeId, 'H5T_NATIVE_FLOAT') == 1);
                    otherwise
                        assert(H5T.equal(typeId, 'H5T_NATIVE_INT') == 1);
                end
                % Get data space
                if ndims(wdata) == 3
                    datadims = fliplr(size(wdata));
                else
                    datadims = fliplr([size(wdata) 1]);
                end
                memspace = H5S.create_simple(3,datadims,[]);
                % Dimensions
                filespace = H5D.get_space(dataset);
                [nDims, h5dims] = H5S.get_simple_extent_dims(filespace);
                H5S.close (filespace);
                % Extend the dataset
                offset = fliplr(offset);
                assert(nDims == numel(offset));
                newdims = offset + datadims;
                dims = max([h5dims; newdims]);
                H5D.set_extent(dataset,dims);
                filespace = H5D.get_space(dataset);
                % Select a hyperslab in extended portion of dataset
                H5S.select_hyperslab(filespace,'H5S_SELECT_SET',...
                    offset,[],[],datadims);
                % Write data
                H5D.write(dataset,'H5ML_DEFAULT', memspace, filespace, ...
                    'H5P_DEFAULT',wdata);
                % Close resources
                H5D.close(dataset);
                H5S.close(memspace);
                H5S.close (filespace);
            else
                writedataset(obj, datasetName, wdata);
            end
        end
        
        %%
        function dims = getdatasetsize(obj, dsetname)
            %% Read dataset dimensions
            % Arguments:
            %   dsetname - full path from the root group to the dataset
            dset = H5D.open(obj.fileID,dsetname);
            space = H5D.get_space(dset);
            [~,dims] = H5S.get_simple_extent_dims(space);
            H5D.close(dset);
        end
        
        %%
        function nbytes = getbytespersample(obj, dsetname)
            %% Read data type
            % Arguments:
            %   dsetname - full path from the root group to the dataset
            dset = H5D.open(obj.fileID,dsetname);
            type = H5D.get_type(dset);
            nbytes = H5T.get_size(type);
            H5T.close(type);
            H5D.close(dset);
        end
        
        function attr = readattribute(obj, objname, attrname)
            %% Read attribute
            % Arguments:
            %   objname - full path from the root group to object the
            %   attribute is related to. Empty if related to file
            %   attrname - attribute name
            if ~isempty(objname)
                oid = H5O.open(obj.fileID,objname);
            else
                oid = obj.fileID;
            end
            attr_id = H5A.open(oid,attrname);
            attr = H5A.read(attr_id);
            H5A.close(attr_id);
            if ~isempty(objname)
                H5O.close(oid);
            end
        end
        
        function createattribute(obj, objname, attrname, value)
            %% Create a new attribute and write its value
            % Arguments:
            %   objname - full path from the root group to object the
            %   attribute is related to. Empty if related to file
            %   attrname - attribute name
            %   value - attribute value
            if ~isempty(objname)
                oid = H5O.open(obj.fileID,objname,'H5P_DEFAULT');
            else
                oid = obj.fileID;
            end
            % Data type
            switch class(value)
                case 'uint16'
                    datatype = H5T.copy('H5T_NATIVE_UINT16');
                case 'uint8'
                    datatype = H5T.copy('H5T_NATIVE_UINT8');
                case 'int16'
                    datatype = H5T.copy('H5T_NATIVE_INT16');
                case 'int8'
                    datatype = H5T.copy('H5T_NATIVE_INT8');
                case 'double'
                    datatype = H5T.copy('H5T_NATIVE_DOUBLE');
                case 'single'
                    datatype = H5T.copy('H5T_NATIVE_FLOAT');
                otherwise
                    datatype = H5T.copy('H5T_NATIVE_INT');
            end
            
            ndim = ndims(value);
            if ndim > 1
                dims = fliplr(size(value));
            else
                dims = 1;
            end
            unlim = H5ML.get_constant_value('H5S_UNLIMITED');
            
            maxdims = unlim * ones(1,ndim);
            dataspace = H5S.create_simple(ndim,dims,maxdims);
            
            attr_id = H5A.create(oid, attrname, datatype, dataspace,...
                'H5P_DEFAULT');
            H5A.write(attr_id, 'H5ML_DEFAULT', value)
            H5A.close(attr_id);
            if ~isempty(objname)
                H5O.close(oid);
            end
        end
        
    end
    
end