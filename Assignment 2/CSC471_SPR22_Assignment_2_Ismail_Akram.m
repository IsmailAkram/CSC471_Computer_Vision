%===================================================
% Computer Vision Programming Assignment 2
% @Zhigang Zhu
% City College of New York
%===================================================

% Ismail Akram 8834

InputImage = 'IDPicture.bmp';
%OutputImage1 = 'IDPicture_bw.bmp';

C1 = imread(InputImage);
[ROWS COLS CHANNELS] = size(C1); % 250 * 250 *3 (RGB)
figure();
subplot(2,2,1), image(C1), title('Original Image');

% % Splitting bands for RGB
% CR1 =uint8(zeros(ROWS, COLS, CHANNELS)); 
% for band = 1 : CHANNELS, 
%     CR1(:,:,band) = (C1(:,:,1));
% end
% 
% CG1 =uint8(zeros(ROWS, COLS, CHANNELS));
% for band = 1 : CHANNELS, 
%     CG1(:,:,band) = (C1(:,:,2));
% end
% 
% CB1 =uint8(zeros(ROWS, COLS, CHANNELS));
% for band = 1 : CHANNELS,
%     CB1(:,:,band) = (C1(:,:,3));
% end
% 
% %I = (0.299 * CR1) + (0.587 * CG1) + (0.114 * CB1); % Gray scale
% I = uint8(round(sum(C1,3)/3)); % Intensity

% Color mapping b/w
MAP = zeros(256, 3);

for i = 1 : 256,  % a comma means pause 
    for band = 1:CHANNELS,
        MAP(i,band) = (i-1)/255; % scaling 0-255 into 0 to 1
    end 
end

% Transform colour image to intensity image
I = (0.299*C1(:,:,1)) + (0.587*C1(:,:,2)) + (0.114*C1(:,:,3));
subplot(2,2,2), image(I), title('Intensity Image');
colormap(MAP); % maps the image to the color map (b/w)

%% Question 1
% Histogram
Intensityhistogram = zeros(256, 1);
for i = 1:250
    for j = 1:250
        Intensityhistogram(I(i,j) + 1) = ...
            Intensityhistogram(I(i,j)+1)+1; % increment by 1
    end
end

subplot(2,2,[3,4]);
bar(Intensityhistogram, 0.8); % histogram
title('Histogram');

% Contrast enhancement
Bright_Enhancement = I + 30;
figure();
subplot(2,2,1), image(I);
title('Original Intensity');
subplot(2,2,2), image(Bright_Enhancement);
title('Point Transforms: Brightness');
colormap(MAP);

Intensityhistogram = zeros(256, 1); % Histogram for Contrast enhancement
for i = 1:250
    for j = 1:250
        Intensityhistogram(Bright_Enhancement(i, j) + 1) = ...
            Intensityhistogram(Bright_Enhancement(i, j) + 1) + 1;
        % At every intensity of IntensityHistogram(0+1), we increment by 1 
    end
end
subplot(2,2,[3,4]), bar(Intensityhistogram, 0.8);
title('Histogram for brightness');

% Linear Stretch
contrast = ((255/225)*(I-20));
figure();
subplot(2,2,1), image(I);
title('Original Intensity');
subplot(2,2,2), image(contrast);
title('Contrast');
colormap(MAP);

Intensityhistogram = zeros(256, 1); % Histogram for Contrast enhancement
for i = 1:250
    for j = 1:250
        Intensityhistogram(contrast(i, j) + 1) = ...
            Intensityhistogram(contrast(i, j) + 1) + 1;
        % At every intensity of IntensityHistogram(0+1), we increment by 1 
    end
end
subplot(2,2,[3,4]), bar(Intensityhistogram, 0.8);
title('Histogram for Contrast');

% Thresholding
Thresholding = I;
for i = 1:250
    for j = 1:250
        if Thresholding(i,j) < 105
            Thresholding(i,j) = 0;
        else 
            Thresholding(i,j) = 255;
        end
    end
end
figure();
subplot(2,2,1), image(I);
title('Original Intensity');
subplot(2,2,2), image(Thresholding);
title('Thresholding');
colormap(MAP);

Intensityhistogram = zeros(256, 1); 
for i = 1:250
    for j = 1:250
        Intensityhistogram(Thresholding(i, j) + 1) = ...
            Intensityhistogram(Thresholding(i, j) + 1) + 1;
        % At every intensity of IntensityHistogram(0+1), we increment by 1 
    end
end
subplot(2,2,[3,4]), bar(Intensityhistogram, 0.8);
title('Histogram for Threshold');

% Equalization
Equalized_I = uint8(zeros(size(I,1), size(I,2)));

Counter = imhist(I);
CumFreq = cumsum(Counter) / 62500; % MxN = 250x250 = 62500
output = round(CumFreq * 255); % 0-255
Equalized_I(:) = output(I(:) + 1); % + 1 since we start at index 1

figure();
subplot(2,2,1), image(I);
title('Original Intensity');
subplot(2,2,2), image(Equalized_I);
title('Equalization');
colormap(MAP);

Intensityhistogram = zeros(256, 1); % Histogram for equalization
for i = 1:250
    for j = 1:250
        Intensityhistogram(Equalized_I(i, j) + 1) = ...
            Intensityhistogram(Equalized_I(i, j) + 1) + 1;
        % At every intensity of IntensityHistogram(0+1), we increment by 1 
    end
end
subplot(2,2,[3,4]), bar(Intensityhistogram, 0.8);
title('Histogram for Equalization');

%% Question 2

% Kernal:
% Ax = [-1, 1]; Y-direction
% Ay = [-1,     X-direction
%        1];
% Can divide by 1 to reduce noise

input_image = double(I);

% Pre-allocating the gradient_image matric with zeros
gradient_image = zeros(size(input_image));

% 1x2 operator mask
gradient_image_1x2 = gradient_image;

Mx_1x2 = [-1, 1]; % 1x2 x-axis
% My_1x2 = [-1; 1]; % 2x1 y-axis

for i = 1:size(input_image, 1) - 2
    for j = 1:size(input_image, 2) - 2
  
        % Gradient approximations
        Gx = sum(sum(Mx_1x2.*input_image(i:i+1, j:j+1)));
        % Gy = sum(sum(My_1x2.*input_image(i:i+1, j:j+1)));
                 
        % Calculate magnitude of vector
        gradient_image_1x2(i+1, j+1) = sqrt(Gx.^2);         
    end
end

% Displaying Gradient Image
% gradient_image_1x2 = normalize(gradient_image_1x2);
gradient_image_1x2 = uint8(gradient_image_1x2)/4; % normalizing using |-1| + |1| + |-1| + |1| = 4

No8 = figure;
imshow(gradient_image_1x2); 
title('Gradient Image (1x2)');

% 2x1 operator mask
gradient_image_2x1 = gradient_image;

% Mx_2x1 = [-1, 1]; % 1x2 x-axis
My_2x1 = [-1; 1]; % 2x1 y-axis

for i = 1:size(input_image, 1) - 2
    for j = 1:size(input_image, 2) - 2
  
        % Gradient approximations
        % Gx = sum(sum(Mx_2x1.*input_image(i:i+1, j:j+1)));
        Gy = sum(sum(My_2x1.*input_image(i:i+1, j:j+1)));
                 
        % Calculate magnitude of vector
        gradient_image_2x1(i+1, j+1) = sqrt(Gy.^2);  
    end
end

% Displaying Gradient Image
% gradient_image_2x1 = normalize(gradient_image_2x1);
gradient_image_2x1 = uint8(gradient_image_2x1)/4;
No9 = figure;
imshow(gradient_image_2x1); 
title('Gradient Image (2x1)');

% Combined 2x2 operator mask
gradient_image_2x2 = gradient_image;

Mx_2x2 = [-1, 1]; % 1x2 x-axis
My_2x2 = [-1; 1]; % 2x1 y-axis

for i = 1:size(input_image, 1) - 2
    for j = 1:size(input_image, 2) - 2
  
        % Gradient approximations
        Gx = sum(sum(Mx_2x2.*input_image(i:i+1, j:j+1)));
        Gy = sum(sum(My_2x2.*input_image(i:i+1, j:j+1)));
                 
        % Calculate magnitude of vector
        gradient_image_2x2(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end

% Displaying Gradient Image
% gradient_image_2x2 = normalize(gradient_image_2x2);
gradient_image_2x2 = uint8(gradient_image_2x2)/4;
No10 = figure;
imshow(gradient_image_2x2); 
title('Gradient Image (2x2)');

% Sobel Operator Mask
gradient_image_S = gradient_image;

Mx_S = [-1 0 1; -2 0 2; -1 0 1];
My_S = [-1 -2 -1; 0 0 0; 1 2 1];

% Edge Detection Process
% When i = 1 and j = 1, then gradient_image pixel  
% position will be gradient_image(2, 2)
% The mask is of 3x3, so we need to traverse 
% to gradient_image(size(input_image, 1) - 2
%, size(input_image, 2) - 2)
% Thus we are not considering the borders.

for i = 1:size(input_image, 1) - 2
    for j = 1:size(input_image, 2) - 2
  
        % Gradient approximations
        Gx = sum(sum(Mx_S.*input_image(i:i+2, j:j+2)));
        Gy = sum(sum(My_S.*input_image(i:i+2, j:j+2)));
                 
        % Calculate magnitude of vector
        gradient_image_S(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);         
    end
end

% Displaying Gradient Image
% gradient_image_S = normalize(gradient_image_S);
gradient_image_S = uint8(gradient_image_S)/8;

No11 = figure;
imshow(gradient_image_S); 
title('Gradient Image (Sobel 3x3)');

% Subtracting 1x2 image gradient from Sobel
Sobel_1x2_difference = gradient_image_S - gradient_image_1x2;
No12 = figure;
imshow(Sobel_1x2_difference);
title('Difference between Sobel & 1x2');

%% Question 3
% % % No14 = figure;
% % % imhist(gradient_image_2x2); NEEDED?
% % % title('threshold-test')

% You may first generate a histogram of each combined gradient map,  
% and only keep certain percentage of pixels
% (e.g.  5% of the pixels with the highest gradient  values) 
% as edge pixels (edgels)

% Get the histogram values position in original picture. 
% Take that position, divide it by total possible pixels. 
% if equal .95, set that value equal to threhsold.
threshold = 0;
position = zeros(256,1);
var = ((250 * 250) * (1 - 0.95)); 
result = 0;
for i = 1:256
     count = IntensityHistogram(i);
     result = result + count;
     if (result >= var)
         threshold = i;
         break;
     end
end

% disp(threshold);

% 1x2 Edge Map
output_image_1x2 = max(gradient_image_1x2, threshold);
output_image_1x2(output_image_1x2 == round(threshold)) = 0;

% Displaying Output Image
%output_image_1x2 = uint8(output_image_1x2);
output_image_1x2 = im2bw(output_image_1x2);
No13 = figure;
imshow(output_image_1x2); 
title('Edge Detected Image (1x2)');

% 2x1 Edge Map
output_image_2x1 = max(gradient_image_2x1, threshold);
output_image_2x1(output_image_2x1 == round(threshold)) = 0;

% Displaying Output Image
output_image_2x1 = im2bw(output_image_2x1);
No14 = figure;
imshow(output_image_2x1); 
title('Edge Detected Image (2x1)');

% 2x2 Edge Map
output_image_2x2 = max(gradient_image_2x2, threshold);
output_image_2x2(output_image_2x2 == round(threshold)) = 0;

% Displaying Output Image
output_image_2x2 = im2bw(output_image_2x2); % im2bw
No15 = figure;
imshow(output_image_2x2); 
title('Edge Detected Image (2x2)');

% Sobel Edge Map
output_image_S = max(gradient_image_S, threshold);
output_image_S(output_image_S == round(threshold)) = 0;

% Displaying Output Image
output_image_S = im2bw(output_image_S); % im2bw
No16 = figure;
imshow(output_image_S); 
title('Edge Detected Image (Sobel 3x3)');

%% Question 4
% Sobel 5x5 Operator Mask
gradient_image_S5x5 = gradient_image;

Mx_S5x5 = [+2 +2 +4 +2 +2; +1 +1 +2 +1 +1; 0 0 0 0 0; -1 -1 -2 -1 -1; -2 -2 -4 -2 -2];
My_S5x5 = [+2 +1 0 -1 -2; +2 +1 0 -1 -2; +4 +2 0 -2 -4; +2 +1 0 -1 -2; +2 +1 0 -1 -2];

% Edge Detection Process
% When i = 1 and j = 1, then gradient_image pixel  
% position will be gradient_image(2, 2)
% The mask is of 3x3, so we need to traverse 
% to gradient_image(size(input_image, 1) - 2
%, size(input_image, 2) - 2)
% Thus we are not considering the borders.

for i = 1:size(input_image, 1) - 4
    for j = 1:size(input_image, 2) - 4
  
        % Gradient approximations
        Gx = sum(sum(Mx_S5x5.*input_image(i:i+4, j:j+4)));
        Gy = sum(sum(My_S5x5.*input_image(i:i+4, j:j+4)));
                 
        % Calculate magnitude of vector
        gradient_image_S5x5(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end

% Displaying Gradient Image
gradient_image_S5x5 = uint8(gradient_image_S5x5); % normalizing
No16 = figure;
imshow(gradient_image_S5x5); 
title('Gradient Image (Sobel 5x5)');

output_image_S5x5 = max(gradient_image_S5x5, threshold);
output_image_S5x5(output_image_S == round(threshold)) = 0;

% Displaying Output Image
output_image_S5x5 = im2bw(output_image_S5x5); % im2bw
No17 = figure;
imshow(output_image_S5x5); 
title('Edge Detected Image (Sobel 5x5)');

% Sobel 7x7 Operator Mask
gradient_image_S7x7 = gradient_image;

Mx_S7x7 = [+3 +3 +3 +6 +3 +3 +3; +2 +2 +2 +4 +2 +2 +2; +1 +1 +1 +2 +1 +1 +1; 0 0 0 0 0 0 0; -1 -1 -1 -2 -1 -1 -1; -2 -2 -2 -4 -2 -2 -2; -3 -3 -3 -6 -3 -3 -3];
My_S7x7 = [+3 +2 +1 0 -1 -2 -3; +3 +2 +1 0 -1 -2 -3; +3 +2 +1 0 -1 -2 -3; +6 +4 +2 0 -2 -4 -6; +3 +2 +1 0 -1 -2 -3; +3 +2 +1 0 -1 -2 -3; +3 +2 +1 0 -1 -2 -3];

% Edge Detection Process
% When i = 1 and j = 1, then gradient_image pixel  
% position will be gradient_image(2, 2)
% The mask is of 3x3, so we need to traverse 
% to gradient_image(size(input_image, 1) - 2
%, size(input_image, 2) - 2)
% Thus we are not considering the borders.

for i = 1:size(input_image, 1) - 6
    for j = 1:size(input_image, 2) - 6
  
        % Gradient approximations
        Gx = sum(sum(Mx_S7x7.*input_image(i:i+6, j:j+6)));
        Gy = sum(sum(My_S7x7.*input_image(i:i+6, j:j+6)));
                 
        % Calculate magnitude of vector
        gradient_image_S7x7(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end

% Displaying Gradient Image
gradient_image_S7x7 = uint8(gradient_image_S7x7); % normalizing
No18 = figure;
imshow(gradient_image_S7x7); 
title('Gradient Image (Sobel 7x7)');

output_image_S7x7 = max(gradient_image_S7x7, threshold);
output_image_S7x7(output_image_S == round(threshold)) = 0;

% Displaying Output Image
output_image_S7x7 = im2bw(output_image_S7x7); % im2bw
No19 = figure;
imshow(output_image_S7x7); 
title('Edge Detected Image (Sobel 7x7)');

%% Question 5
% Converting to double

% CRBand = uint8(CR1);
% RBand_input_image = double(CRBand); % red band

% Red band
RBand_input_image = C1(:, :, 1);
RBand_input_image = double(RBand_input_image);

% Pre-allocating the gradient_image matric with zeros
gradient_image_RB = zeros(size(RBand_input_image));

% Sobel Operator Mask
Mx_S = [-1 0 1; -2 0 2; -1 0 1];
My_S = [-1 -2 -1; 0 0 0; 1 2 1];

% Edge Detection Process
% When i = 1 and j = 1, then gradient_image pixel  
% position will be gradient_image(2, 2)
% The mask is of 3x3, so we need to traverse 
% to gradient_image(size(input_image, 1) - 2
%, size(input_image, 2) - 2)
% Thus we are not considering the borders.

for i = 1:size(RBand_input_image, 1) - 2
    for j = 1:size(RBand_input_image, 2) - 2
  
        % Gradient approximations
        Gx = sum(sum(Mx_S.*RBand_input_image(i:i+2, j:j+2)));
        Gy = sum(sum(My_S.*RBand_input_image(i:i+2, j:j+2)));
                 
        % Calculate magnitude of vector
        gradient_image_RB(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end

% Displaying Gradient Image
gradient_image_RBand = uint8(gradient_image_RB);
No20 = figure;
imshow(gradient_image_RBand); 
title('Gradient Image Red Band (Sobel 3x3)');

output_image_RB = max(gradient_image_RBand, threshold);
output_image_RB(output_image_RB == round(threshold)) = 0;

% Displaying Output Image
output_image_RB = im2bw(output_image_RB); % im2bw
No21 = figure;
imshow(output_image_RB); 
title('Edge Detected Image Red Band (Sobel 3x3)');

% Green band
GBand_input_image = C1(:, :, 1);
GBand_input_image = double(GBand_input_image);

% Pre-allocating the gradient_image matric with zeros
gradient_image_GB = zeros(size(GBand_input_image));

% Sobel Operator Mask
Mx_S = [-1 0 1; -2 0 2; -1 0 1];
My_S = [-1 -2 -1; 0 0 0; 1 2 1];

% Edge Detection Process
% When i = 1 and j = 1, then gradient_image pixel  
% position will be gradient_image(2, 2)
% The mask is of 3x3, so we need to traverse 
% to gradient_image(size(input_image, 1) - 2
%, size(input_image, 2) - 2)
% Thus we are not considering the borders.

for i = 1:size(GBand_input_image, 1) - 2
    for j = 1:size(GBand_input_image, 2) - 2
  
        % Gradient approximations
        Gx = sum(sum(Mx_S.*GBand_input_image(i:i+2, j:j+2)));
        Gy = sum(sum(My_S.*GBand_input_image(i:i+2, j:j+2)));
                 
        % Calculate magnitude of vector
        gradient_image_GB(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end

% Displaying Gradient Image
gradient_image_GBand = uint8(gradient_image_GB);
No22 = figure;
imshow(gradient_image_GBand); 
title('Gradient Image Green Band (Sobel 3x3)');

output_image_GB = max(gradient_image_GBand, threshold);
output_image_GB(output_image_GB == round(threshold)) = 0;

% Displaying Output Image
output_image_GB = im2bw(output_image_GB); % im2bw
No23 = figure;
imshow(output_image_GB); 
title('Edge Detected Image Green Band (Sobel 3x3)');

% Blue band
BBand_input_image = C1(:, :, 1);
BBand_input_image = double(BBand_input_image);

% Pre-allocating the gradient_image matric with zeros
gradient_image_BB = zeros(size(BBand_input_image));

% Sobel Operator Mask
Mx_S = [-1 0 1; -2 0 2; -1 0 1];
My_S = [-1 -2 -1; 0 0 0; 1 2 1];

% Edge Detection Process
% When i = 1 and j = 1, then gradient_image pixel  
% position will be gradient_image(2, 2)
% The mask is of 3x3, so we need to traverse 
% to gradient_image(size(input_image, 1) - 2
%, size(input_image, 2) - 2)
% Thus we are not considering the borders.

for i = 1:size(BBand_input_image, 1) - 2
    for j = 1:size(BBand_input_image, 2) - 2
  
        % Gradient approximations
        Gx = sum(sum(Mx_S.*BBand_input_image(i:i+2, j:j+2)));
        Gy = sum(sum(My_S.*BBand_input_image(i:i+2, j:j+2)));
                 
        % Calculate magnitude of vector
        gradient_image_BB(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end

% Displaying Gradient Image
gradient_image_BBand = uint8(gradient_image_BB);
No24 = figure;
imshow(gradient_image_BBand); 
title('Gradient Image Blue Band (Sobel 3x3)');

output_image_BB = max(gradient_image_BBand, threshold);
output_image_BB(output_image_BB == round(threshold)) = 0;

% Displaying Output Image
output_image_BB = im2bw(output_image_BB); % im2bw
No25 = figure;
imshow(output_image_BB); 
title('Edge Detected Image Blue Band (Sobel 3x3)');

Full_color_band_gradient = (gradient_image_RBand + gradient_image_GBand + gradient_image_BBand)/3;
No27 = figure;
imshow(Full_color_band_gradient); 
title('Gradient Image Full Band (Sobel 3x3)');

Full_color_band_edge = (output_image_RB + output_image_RB + output_image_BB)/3;
No28 = figure;
imshow(Full_color_band_edge); 
title('Edge Detected Image Full Band (Sobel 3x3)');