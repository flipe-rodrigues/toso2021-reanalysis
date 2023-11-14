function b = padarray(varargin)
%PADARRAY Pad array.
%   B = PADARRAY(A,PADSIZE) pads array A with PADSIZE(k) number of zeros
%   along the k-th dimension of A.  PADSIZE should be a vector of
%   nonnegative integers. When A is categorical, the array is padded with
%   <undefined> category
%
%   B = PADARRAY(A,PADSIZE,PADVAL) pads array A with PADVAL. If A is
%   numeric or logical, PADVAL should be a scalar.
%   When A is categorical input, PADVAL can take one of the following values:
%
%       - Valid category in input data specified as string or character array
%       - missing, which corresponds  to <undefined> category (default)
%
%
%   B = PADARRAY(A,PADSIZE,PADVAL,DIRECTION) pads A in the direction
%   specified by DIRECTION. DIRECTION can be one of the following:
%
%       String or character vector values for DIRECTION
%       'pre'         Pads before the first array element along each
%                     dimension .
%       'post'        Pads after the last array element along each
%                     dimension.
%       'both'        Pads before the first array element and after the
%                     last array element along each dimension.
%
%   By default, DIRECTION is 'both'.
%
%   B = PADARRAY(A,PADSIZE,METHOD,DIRECTION) pads array A using the
%   specified METHOD.  METHOD can be one of the following:
%
%       String or character vector values for METHOD
%       'circular'    Pads with circular repetition of elements.
%       'replicate'   Repeats border elements of A.
%       'symmetric'   Pads array with mirror reflections of itself.
%
%   Class Support
%   -------------
%   When padding with a constant value, A can be numeric, categorical or
%   logical. When padding using the 'circular', 'replicate', or 'symmetric'
%   methods, A can be of any class.  B is of the same class as A.
%
%   Example 1:
%   ----------
%  % Add three elements of padding to the beginning of a vector.  The
%  % padding elements contain mirror copies of the array.
%
%       b = padarray([1 2 3 4],3,'symmetric','pre')
%
%   Example 2:
%   ----------
%  % Add three elements of padding to the end of the first dimension of
%  % the array and two elements of padding to the end of the second
%  % dimension.  Use the value of the last array element as the padding
%  % value.
%
%       B = padarray([1 2; 3 4],[3 2],'replicate','post')
%
%   Example 3:
%   ----------
%   % Add three elements of padding to each dimension of a
%   % three-dimensional array.  Each pad element contains the value 0.
%
%       A = [1 2; 3 4];
%       B = [5 6; 7 8];
%       C = cat(3,A,B)
%       D = padarray(C,[3 3],0,'both')
%
%   See also CIRCSHIFT, IMFILTER, MISSING.

%   Copyright 1993-2019 The MathWorks, Inc.

args = matlab.images.internal.stringToChar(varargin);

[a, method, padSize, padVal, direction, catConverter] = ParseInputs(args{:});

b = padarray_algo(a, padSize, method, padVal, direction);

if ~isempty(catConverter)
    b = catConverter.numeric2Categorical(b);
end

%%%
%%% ParseInputs
%%%
function [a, method, padSize, padVal, direction, catConverter] = ParseInputs(varargin)

narginchk(2,4);

% fixed syntax args
a         = varargin{1};
padSize   = varargin{2};

% default values
method    = 'constant';
padVal    = 0;
direction = 'both';
catConverter = [];
isCategoricalInput = false;

if iscategorical(a)
    categoriesIn = categories(a);
    catConverter = images.internal.utils.CategoricalConverter(categoriesIn);
    a = catConverter.categorical2Numeric(a);
    isCategoricalInput = true;
end

validateattributes(padSize, {'double'}, {'real' 'vector' 'nonnan' 'nonnegative' ...
    'integer'}, mfilename, 'PADSIZE', 2);

% Preprocess the padding size
if (numel(padSize) < ndims(a))
    padSize           = padSize(:);
    padSize(ndims(a)) = 0;
end


if nargin > 2
    
    firstStringToProcess = 3;
    
    if ~ischar(varargin{3})
        % Third input must be pad value.
        
        padVal = varargin{3};
        
        if isCategoricalInput
            % For categorical inputs, padVal of missing corresponds to <undefined> no
            % other numeric value is valid
            if ~isnumeric(padVal) && ismissing(padVal)
                padVal = 0;
            else
                error(message('images:padarray:invalidPadvalOrDir',2,padVal))
            end
                
        else
            validateattributes(padVal, {'numeric' 'logical'}, {'scalar'}, mfilename, 'PADVAL', 3);
        end
        firstStringToProcess = 4;
        
    end
    
    % Arg 3 and 4, both can take METHOD and DIRECTION interchangeably. This
    % is NOT a bug, the feature was designed this way. One thing to note is
    % ARG 4 cannot be PADVAL
    for k = firstStringToProcess:nargin
        
        validStrings = {'circular' 'replicate' 'symmetric' 'pre' 'post' 'both'};
        
        if isCategoricalInput && k == 3
            validStrings = [categoriesIn', validStrings]; %#ok<AGROW>
            try
                strVal = validatestring(varargin{k}, validStrings);
            catch 
                error(message('images:padarray:invalidPadvalOrDir',k,varargin{k}))
            end
        else
            strVal = validatestring(varargin{k}, validStrings, mfilename, 'METHOD or DIRECTION', k);
        end
        
        switch strVal
            case {'circular' 'replicate' 'symmetric'}
                method = strVal;
                
            case {'pre' 'post' 'both'}
                direction = strVal;
                
            otherwise
                % categorical pad value 
                % Validate case-sensitivity of input padval as
                % validatestring is case-insensitive. Also use varargin{k}
                % instead of strVal, as validatestring always returns lowercase
                % string
                if ~any(strcmp(varargin{k},categoriesIn))
                    error(message('images:padarray:invalidPadvalOrDir',k,varargin{k}))
                end
                padVal = catConverter.getNumericValue(strVal);
                
        end
    end
end

% Check the input array type
if strcmp(method,'constant') && ~(isnumeric(a) || islogical(a))
    error(message('images:padarray:badTypeForConstantPadding'))
end
