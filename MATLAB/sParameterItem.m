classdef channel
    %CHANNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        thru
        NEXT
        FEXT
    end
    
    properties (Dependent = true)
    end
    
    methods
        % constructor
        function self = channel(spitem)
            self.thru = spitem;
            self.NEXT = {};
            self.FEXT = {};
        end
        
        % identifier
        % getters
        % more methods
        
    end
    
end



%% local functions
function SMM = mixedmode(S,form)
    % DMCM2 - Converts single-ended s-parameters for differential pairs to mixed-mode
    % s-parameters.  See get_spar.m help for input format.
    %
    % William R. Peters 4/12/04
    % 
    % USAGE: [SDD,SCC,SDC,SCD,SMM] = mixedmode(S,form)
    %
    % where
    %
    % for form=1 ***************************************************
    %   S = single ended s-parameters (1 = diff pair 1 positive (near end), 3 is diff pair 1
    %     negative, 2 is diff pair 2 pos.(far end) etc.)
    %   SMM = mixed-mode s-parameters (1 = diff pair 1 differential mode, 2 is diff pair 1
    %       common mode, 3 is diff pair 2 differential mode etc.)
    %
    %   1----2		  
    %   3----4	        1 ==== 2
    %             ===> 
    %   5----6	 	    3 ==== 4
    %   7----8          
    %   ...              ...
    %
    % *************************************************************
    % for form=2 **************************************************
    %   S = single ended s-parameters (1 = diff pair 1 positive, 2 is diff pair 1
    %     negative, 3 is diff pair 2 pos. etc.)
    %   SMM = mixed-mode s-parameters (1 = diff pair 1 differential mode, 2 is diff pair 1
    %       common mode, 3 is diff pair 2 differential mode etc.)
    %
    %   1----N/2+1		  
    %   2----N/2+2	      1 ==== N/4+1
    %               ===>
    %   3----N/2+3	 	  2 ==== N/4+2 
    %   4----N/2+4 
    %   ...	              ...
    % **************************************************************
    %

    N = size(S,1)/2;
    k = size(S,3);

    if form==1
        M0 = [1,0,-1;1,0,1]/sqrt(2);
        M = zeros(2*N);
        for ni = 1:2:N,     
            M(2*ni-1:2*ni,2*ni-1:2*ni+1) = M0;  
            M(2*ni+1:2*ni+2,2*ni:2*ni+2) = M0;  
        end
    elseif form==2
        M0 = [1,-1;1,1]/sqrt(2);
        M = zeros(2*N);
        for ni = 1:N
            M(2*ni-1:2*ni,2*ni-1:2*ni) = M0;
        end
    end
    Mi = inv(M);
    SMM = zeros(size(S));
    for ki = 1:k
        SMM(:,:,ki) = M * squeeze(S(:,:,ki)) * Mi;
    end
    % SDD = SMM(1:2:end,1:2:end,:);
    % SDC = SMM(1:2:end,2:2:end,:);
    % SCD = SMM(2:2:end,1:2:end,:);
    % SCC = SMM(2:2:end,2:2:end,:);
end


