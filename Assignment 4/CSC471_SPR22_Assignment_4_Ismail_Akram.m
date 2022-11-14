%% ===========================================================
% CSC I6716 Computer Vision 
% @ Zhigang Zhu, CCNY
% Homework 4 - programming assignment: 
% Fundamental Matrix and Feature-Based Stereo Matching
% 
% Name: Ismail Akram
% ID: 8834
%
% Note: Please do not delete the commented part of the code.
% I am going to use it to test your program
% =============================================================

%% Read in two images 
imgl = imread('pic410.png'); % changed bpm to png
imgr = imread('pic430.png'); 

% display image pair side by side
[ROWS, COLS, CHANNELS] = size(imgl);
disimg = [imgl imgr];
image(disimg);

% You can change these two numbers, 
% but use the variables; I will use them
% to load my data to test your algorithms

% Total Number of control points
Nc = 12;

% Total Number of test points
Nt = 4;

% After several runs, you may want to save the point matches 
% in files (see below and then load them here, instead of 
% clicking many point matches every time you run the program

% load pl.mat pl;
% load pr.mat pr;

% interface for picking up both the control points and 
% the test points

cnt = 1;
hold;

while(cnt <= Nc+Nt)
    
% size of the rectangle to indicate point locations
dR = 50;
dC = 50;

% pick up a point in the left image and display it with a rectangle....
%%% if you loaded the point matches, comment the point picking up (3 lines)%%%
[X, Y] = ginput(1);
Cl = X(1); Rl = Y(1);
pl(cnt,:) = [Cl Rl 1];

% and draw it 
Cl= pl(cnt,1);  Rl=pl(cnt,2); 
rectangle('Curvature', [0 0], 'Position', [Cl Rl dC dR]);

% and then pick up the correspondence in the right image
%%% if you loaded the point matches, comment the point picking up (three lines)%%%

[X, Y] = ginput(1);
Cr = X(1); Rr = Y(1);
pr(cnt,:) = [Cr-COLS Rr 1];

% draw it
Cr=pr(cnt,1)+COLS; Rr=pr(cnt,2);
rectangle('Curvature', [0 0], 'Position', [Cr Rr dC dR]);
plot(Cr+COLS,Rr,'r*');
drawnow;

cnt = cnt+1;
end

%% Student work (1a) NORMALIZATION: Page 156 of the textbook and Ex 7.6
% --------------------------------------------------------------------
% Normalize the coordinates of the corresponding points so that
% the entries of A are of comparable size
% You do not need to do this, but if you cannot get correct
% result, you may want to use this 

% Normalizing within the function since Matlab is giving Syntax Errors

% END NORMALIZATION %%

%% Student work: (1b) Implement EIGHT_POINT algorithm, page 156
% --------------------------------------------------------------------
[ E_Matrix ] = eightpoint( pts_l, pts_r, w, h) % left and right image points, width, height

% Error checking
[n1, c1] = size(pts_l);
[n2, c2] = size(pts_r);
if((c1 ~= 2) || (c2 ~= 2))
    error('Points are not formated with correct number of coordinates.');
end
if((n1 < 8) || (n2 < 8))
    error('There are not enough points to carry out the operation.');
end

% Arranging data for normalization
p_l = transpose([pts_l(1: 8, :), ones(8, 1)]);
p_r = transpose([pts_r(1: 8, :), ones(8, 1)]);
norm1 = getNormMat2d(p_l);
norm2 = getNormMat2d(p_r);

% Normalization
p_l = norm1 * p_l;
p_r = norm2 * p_r;

p_l = transpose(p_l ./ repmat(p_l(3, :), [3, 1]));
p_r = transpose(p_r ./ repmat(p_r(3, :), [3, 1]));

x1 = p_l(:, 1);
y1 = p_l(:, 2);
x2 = p_r(:, 1);
y2 = p_r(:, 2);

% Generate the A matrix
A = [x1.*x2, y1.*x2, x2, ...
       x1.*y2, y1.*y2, y2, ...
       x1,       y1,     ones(8,1)];

% Singular value decomposition of A
[~, ~, V] = svd(A, 0);

% the estimate of F
F_Matrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];
[U, S, V] = svd(F_Matrix);
F_Matrix = U(:, 1) * S(1,1) * transpose(V(:, 1)) + U(:, 2) * S(2,2) * transpose(V(:, 2));


% Undo the coordinate normalization if you have done normalization
F_Matrix = norm2' * F_Matrix * norm1;

% END of EIGHT_POINT

E_Matrix = h' * F_Matrix * w;

%% Draw the epipolar lines for both the controls points and the test
% points, one by one; the current one (* in left and line in right) is in
% red and the previous ones turn into blue

% Find and extract SURF features
% points1 = detectSURFFeatures(imgl);
% points2 = detectSURFFeatures(imgr);
% [f1,vpts1] = extractFeatures(imgl,points1);
% [f2,vpts2] = extractFeatures(imgr,points2);

% I suppose that your Fundamental matrix is F, a 3x3 matrix

% Student work (1c): Check the accuray of the result by 
% measuring the distance between the estimated epipolar lines and 
% image points not used by the matrix estimation.
% You can insert your code in the following for loop

for cnt=1:1:Nc+Nt,
  an = F_Matrix*pl(cnt,:)';
  x = 0:COLS; 
  y = -(an(1)*x+an(3))/an(2);

  x = x+COLS;
  plot(pl(cnt,1),pl(cnt,2),'r*');
  line(x,y,'Color', 'r');
  [X, Y] = ginput(1); %% the location doesn't matter, press mouse to continue...
  plot(pl(cnt,1),pl(cnt,2),'b*');
  line(x,y,'Color', 'b');

end

% Save the corresponding points for later use... see discussions above
% save pr.mat pr;
% save pl.mat pl;

% Save the F matrix in ascii
save F.txt F_Matrix -ASCII

% Student work (1d): Find epipoles using the EPIPOLES_LOCATION algorithm page. 157
% --------------------------------------------------------------------

% save the eipoles 

save eR.txt eRv -ASCII; 
save eL.txt eRv -ASCII; 

% Student work (2). Feature-based stereo matching
% --------------------------------------------------------------------
%% Try to use the epipolar geometry derived from (1) in searching  
% correspondences along epipolar lines in Question (2). You may use 
% a similar interface  as I did for question (1). You may use the point 
% match searching algorithm in (1) (if you have done so), but this 
% time you need to constrain your search windows along the epipolar lines.