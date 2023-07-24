classdef Dual2
    % Dual2 class (Dual numbers as a pair of real numbers)
    %
    % Definition:
    %   A Dual2 of order n is defined as:
    %
    %      (a+b*epsilon)
    %   
    %   where a and b are real numbers and epsilon != 0 is
	%   such that epsilon*epsilon is equal to zero.
    %
    %   An example of Dual2 is therefore  
    %
    %      (3.6-4.2*epsilon)
    %
    %   Constructor for the Dual2 class:
    %
    %   Dual2(re,du)  % builds a Dual2 number, with real part equal to re
    %                 % and dual part equal to du
    %
    %
    %   % EXAMPLES of usage
    %
    %   disp(' ');
    %   d1 = Dual2(3.6, -4.2);
    %
    %   d2 = Dual2(-1.0, 1.0);     % for convenience, also this constructor
    %                              % is supported
    %
    %   d3 = d1+d2;                % computes the sum
    %   disp(d3);                  % will print (2.6 -3.2ε)
    %   
    %
    % See also Dual2Array
    
    % version 1.0
    
    properties
        re double % real part
        du double % dual part
        
        % DUAL2_FORMAT_FOR_REAL_COEFFICIENTS
    end
    
    methods (Static)
      function out = setgetPrefs(varargin)
         persistent Prefs
         
         %if nargin
         %   Prefs = data;
         %end
         
         if isempty(Prefs)
             Prefs.DISPLAY_FORMAT = 0; % ASCII
             % Prefs.DISPLAY_FORMAT = 1; % INTERMEDIATE
             % Prefs.DISPLAY_FORMAT = 2; % LATEX
             
             Prefs.REAL_COEFF_FORMAT = '.4f';
             %Prefs.REAL_COEFF_FORMAT = '%+g';
             %Prefs.REAL_COEFF_FORMAT = '%+.4f';
             
         end
         
         if nargin == 0
             out = Prefs;
         elseif nargin == 1
             out = Prefs.(varargin{1});
         elseif nargin == 2
             Prefs.(varargin{1}) = varargin{2};
             out = Prefs;
         else
             error('setgetPrefs: too many arguments');
         end
             
         % out = Prefs;
      end
   end
    
    methods
        function obj = Dual2(r,d)
            %Dual2 Constructor 
                                    
            if nargin < 1
                r = 0;
            end
            if nargin < 2
                d = 0;
            end
            obj.re = r;
            obj.du = d;                                                                              
        end
        
        function str = char(obj)
            
            Prefs = obj.setgetPrefs();
            
            switch Prefs.DISPLAY_FORMAT
                case 0 % ASCII
                    str = sprintf(['(%', Prefs.REAL_COEFF_FORMAT, ', %+', Prefs.REAL_COEFF_FORMAT], obj.re, obj.du); % QUI scegliamo 17 cifre della mantissa in base 10, che è coerente con una scelta per la BAN_THRESHOLD 10^-17
                    str = [str,'ε)'];                    
                case 1 % INTERMEDIO
                    
                    str = sprintf(['%', Prefs.REAL_COEFF_FORMAT, ' %+', Prefs.REAL_COEFF_FORMAT], obj.re, obj.du); % QUI scegliamo 17 cifre della mantissa in base 10, che è coerente con una scelta per la BAN_THRESHOLD 10^-17
                    str = [str,'ε'];                    
                    
                 case 2 % LATEX
                    str = sprintf(['%', Prefs.REAL_COEFF_FORMAT, ' %+', Prefs.REAL_COEFF_FORMAT], obj.re, obj.du); % QUI scegliamo 17 cifre della mantissa in base 10, che è coerente con una scelta per la BAN_THRESHOLD 10^-17
                    str = [str,'\epsilon'];
                     
                otherwise
                    error('Unknown Dual2 format');
            end
            
        end
        
        function disp(obj)
            %disp This method displays a Dual2 on the screen
            %   This function is build upon member function 'char'
            disp(char(obj));
        end
        
        function re = getReal(obj)
            %disp This method returns the real part of the Dual2 
            re = obj.re;
        end
        
        function du = getDual(obj)
            %disp This method returns the real part of the Dual2 
            du = obj.du;
        end
        
        function d = exp(obj)
            %disp This method returns the exponential of a dual number
            
            % add test if obj is a real number
            
            % add in setgetPref l'ordine a cui fermarsi
            % nello sviluppo do Taylor dell'esponenziale
            
            d = 1;
            fatt = 1;
            pot_obj = 1;
            for i=1:100
                pot_obj = obj*pot_obj;
                d = pot_obj/fatt + d;
                fatt = fatt*i;
            end
            
        end
        
        % BINARY ARITHMETIC OPERATIONS (+,-,*/)
        
        function d3 = plus(d1,d2) % sum of two Dual2 numbers
            if isa(d2,'double')
                if not(isscalar(d2))
                    error('The sum between a scalar Dual2 and a vector/matrix of double is undefined');
                end
                d2 = Dual2(d2); % type promotion
            end
            d3 = Dual2(d1.re+d2.re, d1.du+d2.du);
        end
        
        function d3 = minus(d1,d2) % subtraction of two Dual2 numbers
            if isa(d2,'double')
                if not(isscalar(d2))
                    error('The subtraction between a scalar Dual2 and a vector/matrix of double is undefined');
                end
                d2 = Dual2(d2); % type promotion
            end
            d3 = Dual2(d1.re-d2.re, d1.du-d2.du);
        end
        
        
        function d3 = mtimes(d1,d2) % multiplication
            % mtimes This function computes the multiplication between two Dual2.
            if isa(d2,'double')
                if not(isscalar(d2))
                    error('The multiplication between a scalar Dual2 and a vector/matrix of double is undefined');
                end
               d2 = Dual2(d2);
            end
            d3 = Dual2(d1.re*d2.re, d1.re*d2.du + d1.du*d2.re);
        end
             
        function d3 = mrdivide(d1,d2) % division
            % mtimes This function computes the division between two Dual2.
            if isa(d2,'double')
                if not(isscalar(d2))
                    error('The division between a scalar Dual2 and a vector/matrix of double is undefined');
                end
               d2 = Dual2(d2);
            end
            if d2.re == 0
                error('You cannot divide because the real part of the denominator is 0');
            end
            d3 = Dual2(d1.re/d2.re, (d1.du*d2.re - d1.re*d2.du)/(d2.re^2));
        end  

    end % END Methods
end % END Dual2