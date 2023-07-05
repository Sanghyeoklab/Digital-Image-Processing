%% Data clear
clear all; close all; clc % Clear everything.
%% Prob 1
Original_image = imread('Lenna.png'); 
% Read Original image "Lenna.png".
Original_image = double(rgb2gray(Original_image)); 
% Original image is color, so it need to convert gray image.
partition = sub_images(Original_image,[256,256],[8,8]); 
% Original image is 256 x 256 image, so divide it small patches.
%% Prob 2 : DCT approximate
DCT_sub_images = zeros(size(partition));
% Initialize transform matrix.
for i=1:size(partition,3)
    [DCT_sub_images(:,:,i),D] = DCT(partition(:,:,i));
    % Apply DCT to sub images. Also make DCT matrix. 
end
%% Prob 3
% Element product matrix for keep only upper left triangular region and set
% zeros for lower right region.



Left_triangular = [1,  1,  1,  1,  1,  1,  1,  1;
                   1,  1,  1,  1,  1,  1,  1,  0;
                   1,  1,  1,  1,  1,  1,  0,  0;
                   1,  1,  1,  1,  1,  0,  0,  0;
                   1,  1,  1,  1,  0,  0,  0,  0;
                   1,  1,  1,  0,  0,  0,  0,  0;
                   1,  1,  0,  0,  0,  0,  0,  0;
                   1,  0,  0,  0,  0,  0,  0,  0];
DCT_compress_images = zeros(size(DCT_sub_images));
% Initialize approximate DCT image.
for i=1:size(partition,3)
    DCT_compress_images(:,:,i) = DCT_sub_images(:,:,i).*Left_triangular;
    % Do element product matrix to remove lower right region.
end
%% Prob 4
DCT_reconstruct_image_slice = zeros(size(partition));
% Initialize reconstruct matrix.
for i=1:size(partition,3)
    
    DCT_reconstruct_image_slice(:,:,i) = D'*DCT_compress_images(:,:,i)*D;
    % Apply IDCT at approximate DCT image.
end
DCT_reconstruct_image = paste(DCT_reconstruct_image_slice,[8,8],[256,256]);
% Now paste all slice as original image.
%% Prob 2' : K - L transform approximate
for i=1:size(partition,3)
    Covariance_sub_images(:,:,i) = Corvariance_matrix(partition(:,:,i));
    % Get Covariance matrix for every image partition.
end
Covariance = mean(Covariance_sub_images,3);
% Apply average in covariance matrix. Then you get total covariance matrix.
%% Prob 3'
[eigenVector,eigenValue] = eig(Covariance);

% Get eigenvalue and eigenvector by "eig" function. You don't need eigvalue
% just for check
[temp,index] = sort(diag(abs(eigenValue)),'descend');
eigenVector = eigenVector(:,index);
% For set eigenvalue sort in descent
eigenValue = eigenValue(index,index);
% You don't need eigenvalue, just for fit eigenvector.
%% Prob 4'
dominants = 32;
eigenVector = eigenVector(:,1:dominants);
% Remain only dominats(32) eigenvectors.
Image_mean = mean(partition,3); 
Image_mean = reshape(Image_mean,numel(Image_mean),1);
% Get mean of image partitions. Consider arrange in line. 


% Compress part
for i=1:size(partition,3)
    Image = reshape(partition(:,:,i),numel(partition(:,:,i)),1);
    % Arrange image partitions in line
    Image_KL_compress(:,:,i) = eigenVector'*(Image-Image_mean);
    % In this case eigenvector is row so transpose eigenvector matrix.
    % Matrix product eigenvector and image partition that remove mean of 
    % image.
end


% Reconstruct part
for i=1:size(partition,3)
    Image_KL_reconstruct(:,:,i) = eigenVector*Image_KL_compress(:,:,i) + Image_mean;
    Image_KL_slice(:,:,i) = reshape(Image_KL_reconstruct(:,:,i),8,8);
end
KL_transform = paste(Image_KL_slice,[8,8],[256,256]);
% Now paste all slice as original image.
%% Prob 2'' : SVD approximate
% g = u*D*v'. gg_T = u*(D^2)*u', g_Tg = v*(D^2)*v'. So you can get u in
% eigenvector of gg_T and v in eigenvector of g_Tg and D by square root of
% eigenvalue in gg_T or g_Tg.
big_eigenvalue = 4;
for i=1:size(partition,3)
    gg_T(:,:,i) =  partition(:,:,i)*partition(:,:,i)';
    g_Tg(:,:,i) =  partition(:,:,i)'*partition(:,:,i);
    
    [temp1,temp2] = eig(gg_T(:,:,i));
    [temp3,temp3] = eig(g_Tg(:,:,i));
    [temp_value,index] = sort(diag(temp2),'descend'); 
    index = index(1:big_eigenvalue);
    u(:,:,i) = temp1(:,index); X(:,:,i) = temp2(index,index); v(:,:,i) = temp3(:,index);
    X(:,:,i) = sqrt(X(:,:,i));
    % Sort by descent in eigenvalue and remain only dominant eigenvalues 
    % and eigenvectors.
end
%% Prob 3''
% However when you consider eigenvector, consider only sign different. So
% you have to check relation between u and v. g_Tu = v*D. So you can get v
% by v = g_Tu/D.
for i=1:size(partition,3)
    v(:,:,i) = partition(:,:,i)'*u(:,:,i)/X(:,:,i);
end
%% Prob 4''
for i=1:size(partition,3)
    SVD_partition(:,:,i) = u(:,:,i)*X(:,:,i)*v(:,:,i)';
    % Reconstruct partitions by u*D*v'
end
SVD_Image = paste(SVD_partition,[8,8],[256,256]);
% Now paste all slice as original image.
%% Prob 5
figure(1);
colormap gray;
% Plot for gray scale.
subplot(3,3,2); imagesc(Original_image); title("Original Image")
% Plot Original image for compare others.
subplot(3,3,4); imagesc(DCT_reconstruct_image); title("DCT approximated Image")
% Plot DCT_reconstruct image.
subplot(3,3,5); imagesc(KL_transform); title("KL transform Image")
% Plot K-L reconstruct image.
subplot(3,3,6); imagesc(SVD_Image); title("SVD Image")
% Plot SVD reconstruct image.


% Plot error map that absolute value of difference between reconstruct and
% original image. MAE : mean of absolute value.
subplot(3,3,7); imagesc(abs(DCT_reconstruct_image-Original_image)); 
title("Error map_{DCT}  MAE : " + string(mean(abs(DCT_reconstruct_image(:)-Original_image(:)))))
subplot(3,3,8); imagesc(abs(KL_transform-Original_image)); 
title("Error map_{KL transform}  MAE : " + string(mean(abs(KL_transform(:)-Original_image(:)))))
subplot(3,3,9); imagesc(abs(SVD_Image-Original_image)); 
title("Error map_{SVD Image}  MAE: " + string(mean(abs(SVD_Image(:)-Original_image(:)))))

%% Prob 6
% Repeat same processes to another image
% Brodatz Texture images: http://www.ux.uis.no/~tranden/brodatz.html
clear all;
Original_image = imread('D95_icon.jpg'); 
Original_image = double(rgb2gray(Original_image)); 
Original_image = imresize(Original_image,[40,40]);
% For contain detail random texture smaller than 8 X 8 size. Subsampling
% 75 X 75 to 40 X 40.
Original_size = size(Original_image); Partition_size = [8,8];
partition = sub_images(Original_image,Original_size,Partition_size); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DCT method   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(partition,3)
    [DCT_sub_images(:,:,i),D] = DCT(partition(:,:,i));
end

Left_triangular = [1,  1,  1,  1,  1,  1,  1,  1;
                   1,  1,  1,  1,  1,  1,  1,  0;
                   1,  1,  1,  1,  1,  1,  0,  0;
                   1,  1,  1,  1,  1,  0,  0,  0;
                   1,  1,  1,  1,  0,  0,  0,  0;
                   1,  1,  1,  0,  0,  0,  0,  0;
                   1,  1,  0,  0,  0,  0,  0,  0;
                   1,  0,  0,  0,  0,  0,  0,  0];
               
for i=1:size(partition,3)
    DCT_compress_images(:,:,i) = DCT_sub_images(:,:,i).*Left_triangular;
end
for i=1:size(partition,3)
    DCT_reconstruct_image_slice(:,:,i) = D'*DCT_compress_images(:,:,i)*D;
end
DCT_reconstruct_image = paste(DCT_reconstruct_image_slice,Partition_size,Original_size);
%%%%%%%%%%%%%%%%%%%%%%%%   K-L transform method   %%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(partition,3)
    Covariance_sub_images(:,:,i) = Corvariance_matrix(partition(:,:,i));
end
Covariance = mean(Covariance_sub_images,3);

dominants = 32;
[eigenVector,eigenValue] = eig(Covariance);
[temp,index] = sort(diag(abs(eigenValue)),'descend');
index = index(1:dominants);
eigenVector = eigenVector(:,index);
eigenValue = eigenValue(index,index);


Image_mean = mean(partition,3); Image_mean = reshape(Image_mean,numel(Image_mean),1);
for i=1:size(partition,3)
    Image = reshape(partition(:,:,i),numel(partition(:,:,i)),1);
    Image_KL_compress(:,:,i) = eigenVector'*(Image-Image_mean);
end


for i=1:size(partition,3)
    Image_KL_reconstruct(:,:,i) = eigenVector*Image_KL_compress(:,:,i) + Image_mean;
    Image_KL_slice(:,:,i) = reshape(Image_KL_reconstruct(:,:,i),Partition_size(1),Partition_size(2));
end
KL_transform = paste(Image_KL_slice,Partition_size,Original_size);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SVD method   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
big_eigenvalue = 4;
for i=1:size(partition,3)
    gg_T(:,:,i) =  partition(:,:,i)*partition(:,:,i)';
    g_Tg(:,:,i) =  partition(:,:,i)'*partition(:,:,i);
    [temp1,temp2] = eig(gg_T(:,:,i));
    [temp3,temp3] = eig(g_Tg(:,:,i));
    [temp_value,index] = sort(diag(temp2),'descend'); 
    index = index(1:big_eigenvalue);
    u(:,:,i) = temp1(:,index); X(:,:,i) = temp2(index,index); v(:,:,i) = temp3(:,index);
    X(:,:,i) = sqrt(X(:,:,i));
end

for i=1:size(partition,3)
    v(:,:,i) = partition(:,:,i)'*u(:,:,i)/X(:,:,i);
end

for i=1:size(partition,3)
    SVD_partition(:,:,i) = u(:,:,i)*X(:,:,i)*v(:,:,i)';
end
SVD_Image = paste(SVD_partition,Partition_size,Original_size);

figure(2);
colormap gray;
subplot(3,3,2); imagesc(Original_image); title("Original Image")
subplot(3,3,4); imagesc(DCT_reconstruct_image); title("DCT approximated Image")
subplot(3,3,5); imagesc(KL_transform); title("KL transform Image")
subplot(3,3,6); imagesc(SVD_Image); title("SVD Image")

subplot(3,3,7); imagesc(abs(DCT_reconstruct_image-Original_image)); 
title("Error map_{DCT}  MAE : " + string(mean(abs(DCT_reconstruct_image(:)-Original_image(:)))))
subplot(3,3,8); imagesc(abs(KL_transform-Original_image)); 
title("Error map_{KL transform}  MAE : " + string(mean(abs(KL_transform(:)-Original_image(:)))))
subplot(3,3,9); imagesc(abs(SVD_Image-Original_image)); 
title("Error map_{SVD Image}  MAE: " + string(mean(abs(SVD_Image(:)-Original_image(:)))))

%% Functions that I made
% Function for make sub images.
function partition = sub_images(input_image,input_size,output_size)
    Patch_x = input_size(1)/output_size(1);     Patch_y = input_size(2)/output_size(2);
    % value for patch number.
    partition = zeros(output_size(1),output_size(2),Patch_x*Patch_y);
    % Make initial form for output.
    for x = 1:Patch_x
        for y = 1:Patch_y
           % Concat each patches in 3rd dimension.
           partition(:,:,(x-1)*Patch_y+y) = ...
           input_image(output_size(1)*(x-1)+1:output_size(1)*x,output_size(2)*(y-1)+1:output_size(2)*y); 
        end
    end
    
end

% Function for paste sub images.
function Image = paste(input_image,input_size,output_size)
    Patch_x = output_size(1)/input_size(1);     Patch_y = output_size(2)/input_size(2);
    % value for patch number.
    Image = zeros(output_size(1),output_size(2));
    % Make initial form for output.
    for x = 1:Patch_x
        for y = 1:Patch_y
            % Connect each patches by 3rd dimension.
           Image(input_size(1)*(x-1)+1:input_size(1)*x,input_size(2)*(y-1)+1:input_size(2)*y)...
           = input_image(:,:,(x-1)*Patch_y+y); 
        end
    end
end

% Function for DCT
function [Output_Image,D] = DCT(Input_Image)
    N = size(Input_Image,1);
    [n,m] = meshgrid(0:N-1);
    D = sqrt(2/N) * cos((2*n+1).*m*pi/(2*N));
    % s(m,n) = sqrt(2/N)*c(m)*cos[(2n+1)*pi*m/2N].
    D(1,:) = D(1,:)/sqrt(2);
    % c(0) = 1/sqrt(2), c(m) = 1 in m is not 0.
    Output_Image = D*Input_Image*D';
    % D is not unitary matrix, so transpose in second D.
end

% Function for Corvariance matrix
function matrix = Corvariance_matrix(Input_image)
%     Consider image is ensemble.
    Image = reshape(Input_image,numel(Input_image),1);
    matrix = Image*Image';
end