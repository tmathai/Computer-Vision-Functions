% function [plCoeff, inliersH1,inliersFullPlane1] = findPlaneEq(P1, Pfull, iterations) 
function [plCoeff, X, Y, ZZ] = findPlaneEq(P1, iterations) 
%
% P1 -- N x 3 points
% Pts1, Pts2 - 3 x N points
%
% output
%       plane eqn coeffs - 1 row vector [d, z, y, x]
%       inlier indicies - 1 binary col vector 

numPts = size(P1,1);

% get the 3D points 
x = P1(:,1);
y = P1(:,2);
z = P1(:,3);

P1new = [P1 ones(length(x),1)];

realdot = @(u, v) u*transpose(v);

maxInliers = [];

eqnCoeffs = [];

inliersStore = []; 

inlierThresh = 0.00009; 



avgDV = [];
pEq = [];
cm = [];

for count=1:iterations
    
    % randomly pick four points
    rindex = randperm(numPts,3);
    % extract those four random points from the point list
    pts = P1new(rindex,:);
    
    % get the mean along the column 
    cmm = mean(pts,1);
    
    % subtract the mean from the three points to center each point 
    ptsNew = bsxfun(@minus, pts, cmm);
    
    % compute SVD on the centered three points 
    [~,~,V] = svd(ptsNew,0); 
    
    % the plane equation coefficients is the third column of the V orthogonal matrix 
    P = V(:,3);
    
    % compute the average distance
    avgDist = mean( (abs(P(1)*x + P(2)*y + P(3)*z + P(4) - dot(P,cmm)) ./ sqrt(sum(P(1:3).^2))) );
    
    avgDV = vertcat(avgDV, avgDist );
    pEq = vertcat(pEq, P'); 
    cm = vertcat(cm, cmm);
    
    plCoeff = pEq;
    
%     % store them into separate values 
%     p1 = pts(1,:);
%     p2 = pts(2,:);
%     p3 = pts(3,:);
%     
%     % get the vector between two point pairs (x1,y1,z1) and (x2,y2,z2)
%     p21 = p2 - p1; 
%     p31 = p3 - p1; 
%     
%     normal = cross(p21, p31);
%     
%     syms x y z
%     point = [x,y,z];
%     planeEq = realdot(point-p1,normal);
%     
%     % get the coefficients of the plane equations
%     % this gives normal and direction 
%     out = vpa(coeffs(planeEq));
%     % coeffs are stored in opposite order -- [d, z, y, x]
%     cd = out(1); 
%     cz = out(2); 
%     cy = out(3); 
%     cx = out(4);
%     
%     % use normal and direction to find the plane equation that best fits
%     % all the points
%     rm = cx*P1(:,1) + cy*P1(:,2) + cz*P1(:,3) + cd;
%     
%     % Inliers had final planeEq values less than 0.0001
%     % values which are above the threshold are outliers that don't lie on plane 
%     % obtained threshold by manually looking at each inlier vs outlier kyp in image,
%     % and figuring out its eqiv 3D point, plugging it into the planeEq, and
%     % checking its value. 
%     inlierIdxs = rm < inlierThresh;
%     
%     % store inlier indicies
%     inliersStore(:,1,count) = inlierIdxs;
%     
%     % count the max number of inliers
%     maxInliers(count) = sum(inlierIdxs); 
%     
%     % store plane eqn coefficients
%     eqnCoeffs(count, :) = out; 
    
end

% get the minimum avg dist from the vector 
[val,ind] = min(avgDV); 
% pull out the plane equation corresponding to the min avg dist
peq = pEq(ind,:);
% pull out the centered mean corresponding to the min avg dist
mv = cm(ind,:); 

N_points = 100; 

xmin = min(P1(:, 1));
xmax = max(P1(:, 1));
ymin = min(P1(:, 2));
ymax = max(P1(:, 2));

grid_spacing_y = (ymax-ymin) / N_points;
grid_spacing_x = (xmax-xmin) / N_points;

normal = peq(1:3);
d = peq(4);

[X, Y] = meshgrid(xmin:grid_spacing_x:xmax, ymin:grid_spacing_y:ymax);
ZZ = (dot(peq,mv) - normal(1)*X - normal(2)*Y - d)./normal(3);


% % compute the plane 
% [X,Y] = meshgrid([min(x):(min(x)+max(x))/50:max(x)],[min(y):(min(y)+max(y))/50:max(y)])
% Zn2 = (dot(peq,mv) - peq(1)*X - peq(2)*Y) ./ peq(3);

% % plot the plane and the original data 
% figure; 
% plot3(x,y,z,'.r'); 
% hold on; 
% plot3(X,Y,Zn2,'.g')
% legend('original data','fitted plane'); 

% % get the index of the maximum inlier count
% [~, k] = min(maxInliers);
% 
% % get plane eqn coefficients with given index 
% plCoeff = eqnCoeffs(k,:);
% 
% % get inliers with given index 
% inliersH1 = inliersStore(:,:,k); 
% 
% rmf = plCoeff(4)*Pfull(:,1) + plCoeff(3)*Pfull(:,2) + plCoeff(2)*Pfull(:,3) + plCoeff(1);
% 
% inlierThresh1 = 0.00003; 
% 
% inliersFullPlane1 = rmf < inlierThresh1; 


end