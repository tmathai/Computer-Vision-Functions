function [H2to1] = computeH_norm(points1, points2)

% input points
% points1, points2 -- 3 x N

% find the rows and columns of each point set
[rows, columns] = size(points1);
[rowsxy, columnsxy] = size(points2);

% calculate the centroid of each point set 
cent = mean(points1(1:2,:)')';
cent2 = mean(points2(1:2,:)')';

% subtract the centroid from the corresponding point set
points1n(1,:) = points1(1,:) - cent(1);
points1n(2,:) = points1(2,:) - cent(2);

points2n(1,:) = points2(1,:) - cent2(1);
points2n(2,:) = points2(2,:) - cent2(2);

% calculate the standard deviation of the points
sd1 = std(points1(1:2,:)');
sd1 = sd1 + (sd1==0);

sd2 = std(points2(1:2,:)');
sd2 = sd2 + (sd2==0);

% to scale by a distance of root 2, divide root 2 by std deviation value
scalingFactor1 = sqrt(2)./sd1;
scalingFactor2 = sqrt(2)./sd2;

% create the transformation matrices
trans1 = [ scalingFactor1(1)    0                 -scalingFactor1(1)*cent(1); ...
           0                scalingFactor1(2)     -scalingFactor1(2)*cent(2); ...
           0                0                     1 ];

 trans2 = [ scalingFactor2(1)    0                   -scalingFactor2(1)*cent2(1); ...
           0                scalingFactor2(2)        -scalingFactor2(2)*cent2(2); ...
           0                0                     1 ];
 
       
% apply transformation matrices to the points   
points1 = trans1 * points1;       
points2 = trans2 * points2;

L = [];
intermZero  = zeros(1,3);
     
for k = 1:columns
    
    % in the for loop, take each point and frame the two equations
    % loop over these two equations and then place them in the right order
    % from Hartley and Zisserman, pg 109 algorithm 
    q1 = points1(:,k);
    q2 = points2(:,k);
    L = [ L;
        q1'*q2(3) intermZero -q1'*q2(1)
        intermZero q1'*q2(3) -q1'*q2(2)
    ];

end
    
    % call svd on the resulting D matrix after the for loop exits
    % the mind eig value corresponds to the eig vector
    [~,S,V] = svd(L, 0); 
    S = diag(S);
    
    htemp = V(:,9);

    % take the index of the eig vector and reshape it fit the requirement
    H = reshape(htemp,3,3)';

    % get the real H 
    H = inv(trans2) * H * trans1;
    
    % to do the image warp, take the inverse of trans2, multiply by H and trans1 
    % this will warp the image like we want it 
    H2to1 = H/H(3,3);
    

end