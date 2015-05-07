function [M2, rotMat, tMat, boolM2] = findM2(F, K1, K2, pts1, pts2)
    
% inputs, F [3x3], K1 [3x4], K2[3x4], pts1[Nx2], pts2[Nx2]
% ouput, M2[3x4], boolM2 - boolean; false - M2 not found; true - M2 found

    % compute the fundamental matrix 
    E = K2' * F * K1;
    
    [U,S,V] = svd(E);
    % take the average of the top two eig values
    m = (S(1,1)+S(2,2))/2;
    % enforce singularity by assigning last eig value to 0 
    E = U*[m,0,0;0,m,0;0,0,0]*V';
    % compute SVD again
    [U,~,V] = svd(E);
    % use W and Z matrices from Martial's slides
    W = [0,-1,0;1,0,0;0,0,1];
    Z = [0,1,0;-1,0,0;0,0,0];
    
    T = U*Z*U';
    t = [T(3,2); T(1,3);T(2,1)];

%     % Account for reflection 
%     % We should return rotation matrices with det(R) == 1
%     if (det(U*W*V')<0)
%         W = -W;
%     end
    
    % two possible rotation matrices
    R1 = U*W*V';
    R2 = U*W'*V';

    % two possible translation vectors
    t1 = t;
    t2 = -t;

    % check determinant for reflection 
    if det(R1) < 0
      R1 = -R1;
    end

    if det(R2) < 0
      R2 = -R2;
    end
    
    % set M1 to be identity [3x4] matrix
    %     M1 = eye(3, 4);
%     M1 = [0.02394,    0.99954,    -0.01856,   41.3607; ...
%           -0.99754,   0.02511,     0.06547,   76.8021; ...
%           0.06591,    0.01695,     0.99768,   902.9826];
    om1R = [0.02394,    0.99954, -0.01856; ...
           -0.99754,   0.02511,  0.06547; ...
           0.06591,    0.01695,  0.99768];
    om1t = [41.3607; 76.8021; 902.9826];

    m1R = om1R';
    m1t = -om1R'*t;
    M1 = [m1R m1t];
    

    
    % get the four different possible M2 matrices 
    % M2 = [R,t], [R,-t], [R',t], [R',-t]    
    M2store = zeros(3,4,4);
    M2store(:,:,1) = [R1,t1];
    M2store(:,:,2) = [R1,t2];
    M2store(:,:,3) = [R2,t1];
    M2store(:,:,4) = [R2,t2];

% 
%     % set M1 as per the parameters provided in the file
%     tRg1 = [-0.01625331773280620100 0.98386957700862299000 -0.17814736905031653000;
%             0.97668439268305030000 -0.02252259937820530100 -0.21349550254417543000;
%             -0.21406407160478280000 -0.17746376518636725000 -0.96056399333613396000];
%     tTg1 = [-0.0283090812583 -0.0366442193256 0.529139415773];
    
%     M1 = [tRg1 tTg1'];
    
    % get the four different possible M2 matrices 
    % M2 = [R,t], [R,-t], [R',t], [R',-t]
%     M2store = zeros(3,4,4);
%     M2store(:,:,1) = [U*W*V',U(:,3)./max(abs(U(:,3)))];
%     M2store(:,:,2) = [U*W*V',-U(:,3)./max(abs(U(:,3)))];
%     M2store(:,:,3) = [U*W'*V',U(:,3)./max(abs(U(:,3)))];
%     M2store(:,:,4) = [U*W'*V',-U(:,3)./max(abs(U(:,3)))];
% %     M2store(:,:,1) = [U*W*V',U(:,3)];
% %     M2store(:,:,2) = [U*W*V',-U(:,3)];
% %     M2store(:,:,3) = [U*W'*V',U(:,3)];
% %     M2store(:,:,4) = [U*W'*V',-U(:,3)];
    
    % M2 exists?
    boolM2 = false; 
    
    rotMat = [];
    tMat = [];
    
    for i=1:4
        % using the four M2s, triangulate all points using them 
        P = triangulatePts(K1*M1, pts1, K2*M2store(:,:,i), pts2);

        % find the M2 for which all the Z coords of the 3D point is > 0
        % this means that the point is in front of both cameras
        if all(P(:,3) > 0)            
            % get the right M2
            M2 = M2store(:,:,i);
            rotMat = M2(:,1:3);     % rotation matrix
            tMat = M2(:,4);         % translation
            boolM2 = true;  % M2 does exist
        end
    end

    
    
    % M2 not found - all 3D points (with Z coords) are not in front of camera 
    % M2 does not exist
    if(boolM2 == false)
        % set M2 to a 3x4 matrix of zeros
        M2 = zeros(3,4); 
        rotMat = M2(:,1:3);     % rotation matrix
        tMat = M2(:,4);         % translation
    end

end
