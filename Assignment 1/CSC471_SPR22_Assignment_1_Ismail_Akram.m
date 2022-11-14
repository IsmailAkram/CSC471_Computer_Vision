%===================================================
% Computer Vision Programming Assignment 1
% @Zhigang Zhu, 2003-2009
% City College of New York
%===================================================

% Ismail Akram 8834

% ---------------- Step 1 ------------------------
% Read in an image, get information
% type help imread for more information

InputImage = 'IDPicture.bmp';
OutputImage1 = 'IDPicture_bw.bmp';

C1 = imread(InputImage);
[ROWS COLS CHANNELS] = size(C1); % 250 * 250 *3 (rgb)

% ---------------- Step 2 ------------------------
% If you want to display the three separate bands
% with the color image in one window, here is 
% what you need to do
% Basically you generate three "color" images
% using the three bands respectively
% and then use [] operator to concatenate the four images
% the orignal color, R band, G band and B band

% First, generate a blank image. Using "uinit8" will 
% give you an image of 8 bits for each pixel in each channel
% Since the Matlab will generate everything as double by default
CR1 =uint8(zeros(ROWS, COLS, CHANNELS)); % blank image

% Note how to put the Red band of the color image C1 into 
% each band of the three-band grayscale image CR1
for band = 1 : CHANNELS, % you're making this section into "1" (R)
    CR1(:,:,band) = (C1(:,:,1));
end

% Do the same thing for G
CG1 =uint8(zeros(ROWS, COLS, CHANNELS));
for band = 1 : CHANNELS, % you're making this section into "2" (G)
    CG1(:,:,band) = (C1(:,:,2));
end

% and for B
CB1 =uint8(zeros(ROWS, COLS, CHANNELS));
for band = 1 : CHANNELS, % you're making this section into "3" (B)
    CB1(:,:,band) = (C1(:,:,3));
end

% Whenever you use figure, you generate a new figure window 
No1 = figure;  % Figure No. 1

% This is what I mean by concatenation using []
disimg = [C1, CR1; CG1, CB1]; % the ; starts a new row (i.e. output as 2x2)

% Then "image" will do the display for you!
image(disimg); % figure 1
title('colour image and respective RGB bands')

% ---------------- Step 3 ------------------------
% Now we can calculate its intensity image from 
% the color image. Don't forget to use "uint8" to 
% covert the double results to unsigned 8-bit integers

I1    = uint8(round(sum(C1,3)/3));

% You can definitely display the black-white (grayscale)
% image directly without turn it into a three-band thing,
% which is a waste of memory space

No2 = figure;  % Figure No. 2
image(I1);
title('Incorrect intensity image (using average of RGB no colormap)')

% If you just stop your program here, you will see a 
% false color image since the system need a colormap to 
% display a 8-bit image  correctly. 
% The above display uses a default color map
% which is not correct. It is beautiful, though

% ---------------- Step 4 ------------------------
% So we need to generate a color map for the grayscale
% I think Matlab should have a function to do this,
% but I am going to do it myself anyway.

% Colormap is a 256 entry table, each index has three entries 
% indicating the three color components of the index

MAP = zeros(256, 3);

% For a gray scale C[i] = (i, i, i)
% But Matlab use color value from 0 to 1 
% so I scale 0-255 into 0-1 (and note 
% that I do not use "unit8" for MAP

for i = 1 : 256,  % a comma means pause 
    for band = 1:CHANNELS,
        MAP(i,band) = (i-1)/255; % scaling 0-255 into 0 to 1
    end 
end

%call colormap to enfore the MAP
colormap(MAP);

% I forgot to mention one thing: the index of Matlab starts from
% 1 instead 0.

% Is it correct this time? Remember the color table is 
% enforced for the current one, which is  the one we 
% just displayed.

% You can test if I am right by try to display the 
% intensity image again:

No3 = figure; % Figure No. 3
image(I1);
title('Corrected intensity image (using average of RGB w/ colormap)')

% See???
% You can actually check the color map using 
% the edit menu of each figure window

% ---------------- Step 5 ------------------------
% Use imwrite save any image
% check out image formats supported by Matlab
% by typing "help imwrite
imwrite(I1, OutputImage1, 'BMP');


% ---------------- Step 6 and ... ------------------------
% Students need to do the rest of the jobs from c to g.
% Write code and comments - turn it in both in hard copies and 
% soft copies (electronically)

% Question 3 
% Generate an intensity image I(x,y) and display it. You should use the equation I = 0.299R + 0.587G + 0.114B

R2 = C1(:, :, 1) * .299;
G2 = C1(:, :, 2) * .587;
B2 = C1(:, :, 3) * .114;

I2 = R2 + G2 + B2;

% I2 = (0.299 * CR1) + (0.587 * CG1) + (0.114 * CB1); % Did not work: it produced a grayscale, 
% $I tried to colormap it but then figure 2 was also affected. I'll just leave this as a comment.

No4 = figure;
image(I2);
title('Intensity image using luminance equation');

% Comparision via algorithm

No5 = figure;
image(abs(double(I1)-(double(I2))));
title('Comparison of I1 and I2');

% Question 4 
K = [4 16 32 64]; % respective 
No6 = figure;

for i = 1:length(K)
    Img_thres = I1;
    intervals = round(linspace(1,256,K(i) + 1)); % +1 since we're using 1-256. To account for 0-255. In MATLAB, our first index is 1.
    for j = 1:length(intervals)-1
        Img_thres(Img_thres > intervals(j) & Img_thres < intervals (j+1)) = intervals(j);
    end
    subplot(2, 2, i);
    imshow(Img_thres);
    title('I1 Quantize');
end

% Question 5
K_c = [2 4]; % respective 
No7 = figure;

for i = 1:length(K_c)
    Img_thres = C1;
    intervals = round(linspace(1,256,K_c(i) + 1)); % +1 since we're using 1-256. To account for 0-255. In MATLAB, our first index is 1.
    for j = 1:length(intervals)-1
        Img_thres(Img_thres > intervals(j) & Img_thres < intervals (j+1)) = intervals(j);
    end
    subplot(1, 2, i);
    imshow(Img_thres);
    title('I1 Quantize');
end

% Question 6
C = ((255/(log(255 + 1)))/3)/100; % I,I' = 255 (maximum value), 3 colour bands, 100 is referring to opacity. So C ~= 0.15...

R3 = C1(:, :, 1);
G3 = C1(:, :, 2);
B3 = C1(:, :, 3);

R4 = C*log(1+(double(R3)));
G4 = C*log(1+(double(G3)));
B4 = C*log(1+(double(B3)));

C1 = cat(3, R4, G4, B4);

No8 = figure;
image(C1);
title('Log Quantization of original three-band colour image C1');