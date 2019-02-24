classdef catilesmap < handle
    
    %%
    properties(SetAccess = private, GetAccess = private)
        tiles_h; % Number of rows of tiles
        tiles_w; % Number of column of tiles
        current_h; % row position of current tile in the map
        current_w; % row position of current tile in the map
    end
    
    %%
    methods(Access = public)
        %%
        function obj = catilesmap(rows,cols)
            %%
            %
            %
            
            %%
            assert(numel(rows) == 1 && numel(cols) == 1);
            obj.tiles_h = rows;
            obj.tiles_w = cols;
            obj.current_h = 1;
            obj.current_w = 1;
        end
        
        %%
        function reinit(obj)
            %%
            %
            %
            
            %%
            obj.current_h = 1;
            obj.current_w = 1;
        end
        
        %%
        function [row,col] = currenttile(obj)
            row = obj.current_h;
            col = obj.current_w;
        end
        
         %%
        function [row,col] = centraltile(obj)
            %%
            %
            %
            
            %%
            row = fix(obj.tiles_h /2);
            col = fix(obj.tiles_w /2);
        end
        
        %%
        function [row,col] = south(obj)
            %%
            %
            %
            
            %%
            if obj.current_h == obj.tiles_h
                row = 0;
                col = 0;
            else
                row = obj.current_h + 1;
                col = obj.current_w;
            end
        end
        
        function [row,col] = north(obj)
            %%
            %
            %
            
            %%
            if obj.current_h == obj.tiles_h
                row = 0;
                col = 0;
            else
                row = obj.current_h + 1;
                col = obj.current_w;
            end
        end
        
        %%
        function [row,col] = east(obj)
            %%
            %
            %
            
            %%
            if obj.current_w == obj.tiles_w
                row = 0;
                col = 0;
            else
                row = obj.current_h;
                col = obj.current_w + 1;
            end
        end
        
        %%
        function [row,col] = west(obj)
            %%
            %
            %
            
            %%
            row = obj.current_h;
            col = obj.current_w - 1;
        end
        
        %%
        function [row, col] = next(obj)
            %%
            %
            %
            
            %%
            row = obj.current_h;
            col = obj.current_w;
            
            if row == 1 || col == obj.tiles_w
                row = col + row;
                col = 1;
                if row > obj.tiles_h
                    e = row - obj.tiles_h;
                    row = row - e;
                    col = col + e;
                end
            else
                row = row - 1;
                col = col + 1;
            end
            
            if row == obj.tiles_h && col == obj.tiles_w
                row = 1;
                col = 1;
            end
            
            obj.current_h = row;
            obj.current_w = col;
        end
        
        %%
        function lindex = twosub2lindex(obj, row, col)
            lindex = (col -1) * obj.tiles_h + row;
        end
        
        %%
        function [row, col] = lindex2twosub(obj, lindex)
            row = rem(lindex, obj.tiles_h);
            col = floor(lindex/obj.tiles_h) + 1;
            if row == 0
                row = obj.tiles_h;
                col = col - 1;
            end
        end
        
        %%
        function nbr = tilesnumber(obj)
            nbr = obj.tiles_h * obj.tiles_w;
        end
        
    end
    
end