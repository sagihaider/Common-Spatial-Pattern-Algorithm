function [result] = f_CSP(varargin)
% This code is for calulating the projection matrix for CSP 
% Haider Raza, Intelligent System Research Center, University of Ulster, Northern Ireland, UK.
%     Raza-H@email.ulster.ac.uk
%     Date: 03-Oct-2014

% Input:
              
%        left:  left hand data
%       right: right hand data
% 
% Output:
%       left: Left hand data
%       right: right hand data    

    if (nargin ~= 2)
        disp('Must have 2 classes for CSP!')
    end
    
    Rsum=0;
    %finding the covariance of each class and composite covariance
    for i = 1:nargin 
        %mean here?
        R{i} = ((varargin{i}*varargin{i}')/trace(varargin{i}*varargin{i}'));%instantiate me before the loop!
        %Ramoser equation (2)
        Rsum=Rsum+R{i};
    end
   
    %   Find Eigenvalues and Eigenvectors of RC
    %   Sort eigenvalues in descending order
    [EVecsum,EValsum] = eig(Rsum);
    [EValsum,ind] = sort(diag(EValsum),'descend');
    EVecsum = EVecsum(:,ind);
    
    %   Find Whitening Transformation Matrix - Ramoser Equation (3)
        W = sqrt(inv(diag(EValsum))) * EVecsum';
    
    
    for k = 1:nargin
        S{k} = W * R{k} * W'; %       Whiten Data Using Whiting Transform - Ramoser Equation (4)
    end
   
    %generalized eigenvectors/values
    [B,D] = eig(S{1},S{2});
    % Simultanous diagonalization
			% Should be equivalent to [B,D]=eig(S{1});
    
    [D,ind]=sort(diag(D));B=B(:,ind);
    
    %Resulting Projection Matrix-these are the spatial filter coefficients
    result = B'*W;
end
