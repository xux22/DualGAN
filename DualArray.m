classdef DualArray
    % DualArray class (Bidimensional Array of Dual2 Numbers)
    %
    %   Constructors for the DualArray class:
    %
    %   DualArray(mat)  % builds a DualArray having the same number of rows
    %                   % and columns of mat, where the element (i,j) is
    %                   % exactly the value of mat(i,j), but represented as
    %                   % a Dual2 of the predetermined degree.
    %
    %   % EXAMPLES of usage
    %
    %   Dual2.setgetPrefs('DISPLAY_FORMAT',0); % ASCII
    %   disp(' ');
    %   mat1 = [ 1 2; 4 5; 6,7]
    %   da1  = DualArray(mat1);
    %   disp(da1); fprintf('\n');
    %
    %   da2 = zeros(2, 3, 'like', DualArray); % this is a 2-by-3 matrix of Dual2 numbers,
    %                                         % all equal to zero.
    %   disp(da2); fprintf('\n');
    %
    %   da3 = eye(3, 'like', DualArray);      % this is a 3-by-3 identity matrix of Dual2
    %   disp(da3); fprintf('\n');
    %
    %
    %   da4 = randn(4, 1, 'like', DualArray); % this is a column vector of random Dual2,
    %                                         % following a normal distribution
    %   disp(da4); fprintf('\n');
    %
    %   da5 = DualArray(randn(5,1));              % this is a 5-by-1 matrix
    %                                             % of Dual2 numbers, with
    %                                             % real part random and dual
    %                                             % part equals to 0
    %                                        
    %   disp(da5); fprintf('\n');
    %
    %   da6 = DualArray(randn(5,3), ones(5,3));   % this is a 5-by-3 matrix of Dual numbers 
    %                                             % with Real part random 
    %                                             % and Dual part equals to
    %                                             % 1
    %                                            
    %   disp(da6); fprintf('\n');
    %
    %   co7 = complex(randn(5,3), ones(5,3));      % you can build a Dual Matrix
    %                                              % with a Matrix of
    %                                              % Complex numbers
    %   da7 = DualArray(co7);                 
    %                                         
    %   disp(da7); fprintf('\n');    
    %
    % See also Dual2, Dual2.setgetPrefs(...)
    
    properties
        dArr % A DualArray is a bidimensional array of Dual2
    end
    methods (Hidden)
        % ------- Redefine the zeros -------
        function obj = zerosLike(obj,varargin)
            if nargin == 1
                error('Please provide the size');
            end
            % With 1-dim, considered a Squared Matrix
            if nargin == 2
                dim = varargin{1};
                obj(dim,dim) = obj;
                for r = 1:dim
                    for c = 1:dim
                        obj(r,c).dArr = Dual2(0);
                    end
                end
            end
            % With 2-dim, the user define both dimensions
            if nargin == 3
                nR = varargin{1};
                nC = varargin{2};
                obj(nR,nC) = obj;
                for r = 1:nR
                    for c = 1:nC
                        obj(r,c).dArr = Dual2(0);
                    end
                end
            end
        end % END zeros
        % ------- Redefine the ones -------
        function obj = onesLike(obj,varargin)
            if nargin == 1
                error('Please provide the size');
            end
            % With 1-dim, considered a Squared Matrix
            if nargin == 2
                dim = varargin{1};
                obj(dim,dim) = obj;
                for r = 1:dim
                    for c = 1:dim
                        obj(r,c).dArr = Dual2(1);
                    end
                end
            end
            % With 2-dim, the user define both dimensions
            if nargin == 3
                nR = varargin{1};
                nC = varargin{2};
                obj(nR,nC) = obj;
                for r = 1:nR
                    for c = 1:nC
                        obj(r,c).dArr = Dual2(1);
                    end
                end
            end
        end % END ones
        % ------- Redefine the eye -------
        function obj = eyeLike(obj,varargin)
            if nargin == 1
                error('Please provide the size');
            end
            if nargin == 2 || ( nargin == 3 && ( varargin{1} == varargin{2}) )
                dim = varargin{1};
                obj(dim,dim) = obj;
                for r = 1:dim
                    for c = 1:dim
                        if r == c
                            obj(r,c).dArr = Dual2(1);
                        else
                            obj(r,c).dArr = Dual2(0);
                        end
                    end
                end
            else
                error('Eye requires that both dimensions are equal');
            end
        end % END eye
        % ------- Redefine the rand -------
        function obj = randLike(obj,varargin)
            if nargin == 1
                error('Please provide the size');
            end
            
            if nargin == 2
                dim = varargin{1};
                obj(dim,dim) = obj;
                for r = 1:dim
                    for c = 1:dim
                        obj(r,c).dArr = Dual2(rand(1,1));
                    end
                end
            end
            if nargin == 3
                nR = varargin{1};
                nC = varargin{2};
                obj(nR,nC) = obj;
                for r = 1:nR
                    for c = 1:nC
                        obj(r,c).dArr = Dual2(rand(1,1));
                    end
                end
            end
        end % END rand
        % ------- Redefine the randn -------
        function obj = randnLike(obj,varargin)
            if nargin == 1
                error('Please provide the size');
            end

            if nargin == 2
                dim = varargin{1};
                obj(dim,dim) = obj;
                for r = 1:dim
                    for c = 1:dim  
                        obj(r,c).dArr = Dual2(randn(1,1));
                    end
                end
            end
            if nargin == 3
                nR = varargin{1};
                nC = varargin{2};
                obj(nR,nC) = obj;
                for r = 1:nR
                    for c = 1:nC
                        obj(r,c).dArr = Dual2(randn(1,1));
                    end
                end
            end
        end % END randn
    end
    
    methods
        % Constructor
        function obj = DualArray(mat, mat_du)
            % DualArray constructor for DualArray
            if nargin ~= 0
                if nargin < 2
                    mat_du = zeros(size(mat));
                end
                if any( size(mat) ~= size(mat_du))
                    error('The matrices for the real part and dual part must have same size');
                end
                [nR,nC] = size(mat);
                obj(nR,nC) = obj;
                %obj = repmat(obj, nR, nC);
                for r = 1:nR
                    for c = 1:nC
                        if isreal(mat)
                           obj(r,c).dArr = Dual2(mat(r,c), mat_du(r,c)); 
                        else % assuming complex
                           obj(r,c).dArr = Dual2(real(mat(r,c)), imag(mat(r,c)));
                        end
                    end
                end
            end
        end % constructor
        
        function str = char(obj)
            error('To be implemented...');
        end % char
        
        function disp(obj)
            [nR, nC] = size(obj);
            for r = 1:nR
                for c = 1:nC
                    fprintf('%s',char(obj(r,c).dArr));
                    if c~=nC
                        fprintf(' , ');
                    end
                end
                fprintf('\n');
            end
        end % disp
        
%         function c = double(obj)
%             R = size(obj,1);
%             aux = obj(1,1).dArr;
%             c = aux.coef;
%             for i=2:R
%                 aux = obj(i,1).dArr;
%                 c = [c; aux.coef];
%             end
%         end % double

        function d2 = getAsDual2(obj, row, column)
             d2 = obj(row,column).dArr;
        end     
        
        function mat_re = getReal(obj)
            mat_re = zeros(size(obj));
            
            nR = size(obj,1);
            nC = size(obj,2);
            
            for r=1:nR
                for c=1:nC
                    d2 = getAsDual2(obj,r,c);
                    mat_re(r,c) = getReal(d2);                    
                end
            end
        end

        function mat_du = getDual(obj)
            mat_du = zeros(size(obj));
            
            nR = size(obj,1);
            nC = size(obj,2);
            
            for r=1:nR
                for c=1:nC
                    d2 = getAsDual2(obj,r,c);
                    mat_du(r,c) = getDual(d2);                    
                end
            end
        end
 
        function abs_as_a_new_DualDarray = abs(obj)
            vec_re = zeros(size(obj));
            vec_du = zeros(size(obj));
            
         
            for i=1:length(obj)
                vec_re(i) = abs(getReal(obj(i)));
                %vec_du(i) = getDual(obj(i));
                vec_du(i) = sign(getReal(obj(i)));
            end
            abs_as_a_new_DualDarray = DualArray(vec_re, vec_du);
        end
 
        function b = subsref(obj,s)
            switch s(1).type
                case '()'
                    ind = s.subs{:};
                    if length(ind) > 1
                        b = obj(ind,1); % returning a DualArray
                    else
                        b = obj(ind,1).dArr; % returning a Dual2
                    end
                 case '{}'
                    ind = s.subs{:};
                    if length(ind) > 1
                        error ('not supported');
                    else
                        b = obj(ind,1).dArr; % returning a Dual2
                    end
                otherwise
                    error('Specify value for x as obj(x)')
            end
        end % subsref
        
        %-BEGIN ARITHMETIC OPERATIONS ----------------------------------------
        % Aritmetic operations between two DualArray or between
        % a DualArray and a double/Dual2
        
        function dArr3 = plus(dArr1, dArr2) % element-wise addition
            if ( size(dArr1,2) > 1 ) % check if the first argument is not a matrix
                error('This function is only available for column vectors');
            end
            dArr3 = zeros(size(dArr1,1), 1, 'like',DualArray);
            if isscalar(dArr2)
                if isa(dArr2,'double')
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr + dArr2; % addition of a Dual2 by a double
                    end
                else
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr + dArr2;
                    end
                end
            else
                if size(dArr1,1) ~= size(dArr2,1)
                    error ('The two DualArray must have the same length');
                else
                    for r=1:size(dArr1,1)
                        if isa(dArr2(r),'double')
                            dArr3(r).dArr = dArr1(r).dArr + dArr2(r);
                        else
                            dArr3(r).dArr = dArr1(r).dArr + dArr2(r).dArr;
                        end
                    end
                end
            end
        end % plus
        
        function dArr3 = minus(dArr1, dArr2) % element-wise subtraction
            if ( size(dArr1,2) > 1 ) % check if the first argument is not a matrix
                error('This function is only available for column vectors');
            end
            dArr3 = zeros(size(dArr1,1), 1, 'like',DualArray);
            if isscalar(dArr2)
                if isa(dArr2,'double')
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr - dArr2; % subtraction of a Dual2 by a double
                    end
                else
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr - dArr2;
                    end
                end
            else
                if size(dArr1,1) ~= size(dArr2,1)
                    error ('The two DualArray must have the same length');
                else
                    for r=1:size(dArr1,1)
                        if isa(dArr2(r),'double')
                            dArr3(r).dArr = dArr1(r).dArr - dArr2(r);
                        else
                            dArr3(r).dArr = dArr1(r).dArr - dArr2(r).dArr;
                        end
                    end
                end
            end
        end % minus
        
        function dArr3 = times(dArr1, dArr2) % element-wise product
            
            if ( size(dArr1,2) > 1 ) % check if the first argument is not a matrix
                error('This function is only available for column vectors');
            end
            dArr3 = zeros(size(dArr1,1), 1, 'like',DualArray);
            if isscalar(dArr2)
                if isa(dArr2,'double')
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr * dArr2; % multiplication of a Dual2 by a double
                    end
                else
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr * dArr2;
                    end
                end
            else
                if size(dArr1,1) ~= size(dArr2,1)
                    error ('The two DualArray must have the same length');
                else
                    for r=1:size(dArr1,1)
                        if isa(dArr2(r),'double')
                            dArr3(r).dArr = dArr1(r).dArr * dArr2(r);
                        else
                            dArr3(r).dArr = dArr1(r).dArr * dArr2(r).dArr;
                        end
                    end
                end
            end
        end % mtimes
        
        function dArr3 = rdivide(dArr1, dArr2) % element-wise division
            if ( size(dArr1,2) > 1 ) % check if the first argument is not a matrix
                error('This function is only available for column vectors');
            end
            dArr3 = zeros(size(dArr1,1), 1, 'like',DualArray);
            if isscalar(dArr2)
                if isa(dArr2,'double')
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr / dArr2; % division of a Dual2 by a double
                    end
                else
                    for r=1:size(dArr1,1)
                        dArr3(r).dArr = dArr1(r).dArr / dArr2;
                    end
                end
            else
                if size(dArr1,1) ~= size(dArr2,1)
                    error ('The two DualArray must have the same length');
                else
                    for r=1:size(dArr1,1)
                        if isa(dArr2(r),'double')
                            dArr3(r).dArr = dArr1(r).dArr / dArr2(r);
                        else
                            dArr3(r).dArr = dArr1(r).dArr / dArr2(r).dArr;
                        end
                    end
                end
            end
        end % mrdivide
        
        %-END ARITHMETIC OPERATIONS ----------------------------------------
        
        %-BEGIN STATISTICAL OPERATIONS -------------------------------------
        function avg = mean(dArr)
            if size(dArr,2) > 1
                error('To Be Implemented!');
            end
            sum = dArr(1).dArr;
            for i = 2:size(dArr,1)
                sum = sum + dArr(i).dArr;
            end
            avg = sum/length(dArr);
        end % mean
        
        function v = var(dArr)
            if size(dArr,2) > 1
                error('To Be Implemented!');
            end
            avg = mean(dArr);
            slacks = dArr-avg;
            squared_slacks = slacks*slacks;
            v = mean(squared_slacks);
        end % var
        
        
        function exp_dArr = exp(dArr1)
            if size(dArr1,2) > 1
                error('To Be Implemented!');
            end
            exp_dArr = dArr1;
            for i = 1:size(dArr1,1)
                exp_dArr(i).dArr = exp(dArr1(i).dArr);                
            end
        end % exp
        
        
    end % methods
end %classdef

