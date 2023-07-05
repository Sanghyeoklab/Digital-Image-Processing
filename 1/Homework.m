%% Data clear
clear all; close all; clc % Clear everything
%% Prob 1
Original_image = imread('Lenna.png'); 
% Read Original image "Lenna.png"
Original_image = double(rgb2gray(Original_image)); 
% Original image is color, so it need to convert gray image
partition = sub_images(Original_image,[256,256],[8,8]); 
% Original image is 256 x 256 image, so divide it small patches
%% Prob 2
for i=1:size(partition,3) 
    % Do DFT in every original image patches and get 8 x 8 DFT unitary matrix
    [F8s(:,:,i),D8] = DFT(partition(:,:,i));
end
%% Prob 3
% S is satisfied matrix that f4 = Sf8(S^T). 
% f4 is subsampled image that 4 x 4 matrix, 
% f8 is original image that 8 x 8 matrix.
S = [0.5,  0.5,  0,    0,    0,    0,    0,    0;...
     0,    0,    0.5,  0.5,  0,    0,    0,    0;...
     0,    0,    0,    0,    0.5,  0.5,  0,    0;...
     0,    0,    0,    0,    0,    0,    0.5,  0.5];
% Get subsampled image in every patches. 
% Do DFT in every subsampled image. Also get 4 x 4 DFT unitary matrix. 
for i=1:size(partition,3)
    Subsampling_image(:,:,i) = S*partition(:,:,i)*S';
    [F4s(:,:,i),D4] = DFT(Subsampling_image(:,:,i));
end
%% Prob 4
% ^H is hermitian matrix that conjugate after transpose original matrix
% Get subsample matrix in frequency version Sf = D4*S*(D8^H).
S_f = D4*S*D8'; 
for i=1:size(partition,3)
    % Apply subsample matrix in frequency version in every original patches
    Subsample_frequency_version(:,:,i) = S_f*F8s(:,:,i)*S_f.';
end
%% Prob 5
% I is satisfied matrix that f8 = Sf4(Sf^T). 
% f4 is original image that 4 x 4 matrix, 
% f8 is interpolation image that 8 x 8 matrix.
I = [1,      0,      0,      0;...
     1/2,    1/2,    0,      0;...
     0,      1,      0,      0;...
     0,      1/2,    1/2,    0;...
     0,      0,      1,      0;...
     0,      0,      1/2,    1/2;...
     0,      0,      0,      1;...
     0       0       0       1]; 
% Get interpolation image in every patches. 
for i=1:size(partition,3)
    Interpolation_image(:,:,i) = I*Subsampling_image(:,:,i)*I';
end
%% Prob 6
% Get interpolation matrix in frequency version If = D8*I*(D4^H).
I_f = D8*I*D4'; 
for i=1:size(partition,3)
    Interpolate_frequency_version(:,:,i) = I_f*F4s(:,:,i)*I_f.';
end
%% Prob 7
% Apply inverse DFT all matrix that get from Prob4 and Prob 6
for i=1:size(partition,3)
   Subsample_image_by_frequency(:,:,i) = D4'* Subsample_frequency_version(:,:,i) * D4';
   Interpolation_image_by_frequency(:,:,i) = D8'* Interpolate_frequency_version(:,:,i) * D8';
end
% Connect each patches for show image
Subsampling_image = paste(Subsampling_image,[4,4],[128,128]);
Interpolation_image = paste(Interpolation_image,[8,8],[256,256]);
Subsample_image_by_frequency = paste(Subsample_image_by_frequency,[4,4],[128,128]);
Interpolation_image_by_frequency = paste(Interpolation_image_by_frequency,[8,8],[256,256]);

% Consider subsample image is 128 x 128 and others are 256 x 256, fit size
% by fill zeros in background.
% Subsampling_image(129:256,:) = 0;
% Subsampling_image(:,129:256) = 0;
% Subsample_image_by_frequency(129:256,:) = 0;
% Subsample_image_by_frequency(:,129:256) = 0;


% % Show each image
% figure(1)
% colormap gray; % For show gray image
% subplot(2,2,1); imagesc(real(Subsampling_image));
% title("Subsample Image");
% subplot(2,2,2); imagesc(real(Interpolation_image));
% title("Interpolation Image");
% subplot(2,2,3); imagesc(real(Subsample_image_by_frequency));
% title("Frequency Subsample Image");
% subplot(2,2,4); imagesc(real(Interpolation_image));
% title("Frequency Interpolation Image");



% Show each image
Boundary = 0.1; size = 0.18;
figure(1)
colormap gray; % For show gray image
subplot(2,2,1); 
set(subplot(2,2,1),'position',[Boundary,1-Boundary-size,size,size]);
imagesc(real(Subsampling_image));
title("Subsample Image");
subplot(2,2,2); 
set(subplot(2,2,2),'position',[0.45+Boundary,1-Boundary-size*2,size*2,size*2]);
imagesc(real(Interpolation_image));
title("Interpolation Image");
subplot(2,2,3); 
set(subplot(2,2,3),'position',[Boundary,0.45-size,size,size]);
imagesc(real(Subsample_image_by_frequency));
title("Frequency Subsample Image");
subplot(2,2,4); 
set(subplot(2,2,4),'position',[0.45+Boundary,0.45-size*2,size*2,size*2]);
imagesc(real(Interpolation_image));
title("Frequency Interpolation Image");


%% Error between apply image domain and frequency domain
fprintf("Maximum error intensity in subsample method: %f\n\n",...
         max(abs(Subsampling_image(:)-Subsample_image_by_frequency(:))))
fprintf("Maximum error intensity in interpolation method : %f\n\n",...
         max(abs(Interpolation_image(:)-Interpolation_image(:))))


%% Functions that I made
% Function for make sub images
function partition = sub_images(input_image,input_size,output_size)
    Patch_x = input_size(1)/output_size(1);     Patch_y = input_size(2)/output_size(2);
    % value for patch number
    partition = zeros(output_size(1),output_size(2),Patch_x*Patch_y);
    % Make initial form for output
    for x = 1:Patch_x
        for y = 1:Patch_y
           % Concat each patches in 3rd dimension
           partition(:,:,(x-1)*Patch_y+y) = ...
           input_image(output_size(1)*(x-1)+1:output_size(1)*x,output_size(2)*(y-1)+1:output_size(2)*y); 
        end
    end
    
end

% Function for paste sub images
function Image = paste(input_image,input_size,output_size)
    Patch_x = output_size(1)/input_size(1);     Patch_y = output_size(2)/input_size(2);
    % value for patch number
    Image = zeros(output_size(1),output_size(2));
    % Make initial form for output
    for x = 1:Patch_x
        for y = 1:Patch_y
            % Connect each patches by 3rd dimension
           Image(input_size(1)*(x-1)+1:input_size(1)*x,input_size(2)*(y-1)+1:input_size(2)*y)...
           = input_image(:,:,(x-1)*Patch_y+y); 
        end
    end
end

% Function for DFT
function [Output_Image,D] = DFT(Input_Image)
    % Check number of dimension. Consider Input_image is square matrix
    N = size(Input_Image,1);
    D = zeros(N,N);
    for k=1:N
        D(k,:) = (k-1)*(0:1:N-1)/N;
    end
    D = exp(-1j*2*pi*D)/sqrt(N);
    % Consider DFT matrix is unitary matrix. So divide by sqrt(N)
    Output_Image = D*Input_Image*D;
end