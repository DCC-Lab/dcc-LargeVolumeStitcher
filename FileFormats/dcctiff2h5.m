
function dcctiff2h5
%%
[filenames, path] = uigetfile('*.tif', 'Pick files', 'Multiselect','on');

name = filenames{1};
delimiters = strfind(name, '-');
prename = name(1:delimiters(end-2));
postname = name(delimiters(end) + 4:end);

h5name = [path 'h5file.h5'];
h5id = H5F.create(h5name,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

gname = sprintf('/R%02d', 0);
gid= H5G.create(h5id,gname,'H5P_DEFAULT','H5P_DEFAULT',...
    'H5P_DEFAULT');

depth = 1;
row   = 1;
col   = 1;

tic
while(1)
    while(1)
        while(1)
            filename = [prename ...
                sprintf('%03d-%03d-%03d', row-1, col-1, depth-1) ...
                postname];
            if ~isempty(find(strcmp(filenames, filename),1))
                data(:,:,depth) = imread([path filename]);
                depth = depth + 1;
            else
                break;
            end
        end
        
        if depth ~= 1
            sgname = sprintf('C%02d', col-1);
            sgid= H5G.create(gid,sgname,'H5P_DEFAULT','H5P_DEFAULT',...
            'H5P_DEFAULT');
            dims = fliplr(size(data));
            dataspace = H5S.create_simple(length(dims), dims, []);
            dname = 'D0000';
            datatype = H5T.copy('H5T_NATIVE_UINT8');
            did = H5D.create(sgid, dname, datatype,dataspace,'H5P_DEFAULT');
            H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',...
                'H5P_DEFAULT',data);
            H5D.close(did);
            H5G.close(sgid);
            col = col + 1;
            depth = 1;
        else
            break;
        end
    end
    
    H5G.close(gid);
    if col ~= 1
        row = row + 1;
        gname = sprintf('/R%02d', row-1);
        gid= H5G.create(h5id,gname,'H5P_DEFAULT','H5P_DEFAULT',...
            'H5P_DEFAULT');
        col = 1;
    else
        break;
    end
end

H5F.close(h5id);
toc

end