function output = cartesian(varargin)
%   CARTESIAN takes the cartesian product (combinations of) the input
%
%   output = CARTESIAN(A) computes the cartesian product of the rows in A
%
%   output = CARTESIAN(a, b, c, ...) computes the cartesian product of the
%   vectors a, b, c ...
%   
%   INPUT: 
%           varargin: Vectors for which to compute cartesian product [Either a M
%           X N matrix, with M the different vectors and N the length of each
%           vector or M vectors (in which case they can be of different length]
%
%   OUTPUT: 
%           output: Contains all possible combinations of the different elements
%           [K x M matrix]
    if nargin == 1
        if size(varargin{1},1) == 1
            % Actually dealing with just a vector, don't actually need to
            % produce the cartesian product. But it is the wrong way around,
            % need to transpose
            output = transpose(varargin{1});
            return
        elseif size(varargin{1},2) == 1
            % Again just dealing with vector, but now can return it instantly
            output = varargin{1};
            return
        else
            % A matrix is given, need to extract the individual vectors
            Nvecs = size(varargin{1},1);

            varargin = mat2cell(varargin{1}',size(varargin{1},2),ones(1,size(varargin{1},1)));  
        end
    elseif nargin > 1
        Nvecs = nargin;
    end
    
    [output{1:Nvecs}] = ndgrid(varargin{:});
    output = reshape(cat(Nvecs+1,output{:}),[],Nvecs);
end
