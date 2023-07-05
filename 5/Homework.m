% HDR : High Dynamic Range
%% Clear everything
clear all; close all; clc;

%% Image load & show
image01 = imread("image01.bmp");
image02 = imread("image02.bmp");
figure(1)
subplot(1,2,1); imagesc(image01); title("Image 1"); axis off;
subplot(1,2,2); imagesc(image02); title("Image 2"); axis off;
%% Prob 1 : Gamma Correction.
image01_gamma1 = Gamma_correction(image01,0.45);
image02_gamma1  = Gamma_correction(image02,0.45);

image01_gamma2 = Gamma_correction(image01,2.2);
image02_gamma2  = Gamma_correction(image02,2.2);
figure(2)
subplot(2,3,1); imagesc(image01); title("Image 1 (Original)"); axis off;
subplot(2,3,2); imagesc(image01_gamma1); title("Image 1 (gamma = 0.45)"); axis off;
subplot(2,3,3); imagesc(image01_gamma2); title("Image 1 (gamma = 2.2)"); axis off;
subplot(2,3,4); imagesc(image02); title("Image 2 (Original)"); axis off;
subplot(2,3,5); imagesc(image02_gamma1); title("Image 2 (gamma = 0.45)"); axis off;
subplot(2,3,6); imagesc(image02_gamma2); title("Image 2 (gamma = 2.2)"); axis off;
%% Prob 2 : Histogram equalization
image01_histogram_equal = Histogram_Equalization(image01);
image02_histogram_equal = Histogram_Equalization(image02);
figure(3)
subplot(2,2,1); imagesc(image01); title("Image 1 (Original)"); axis off;
subplot(2,2,2); imagesc(image01_histogram_equal); title("Image 1 (Histogram equal)"); axis off;
subplot(2,2,3); imagesc(image02); title("Image 2 (Original)"); axis off;
subplot(2,2,4); imagesc(image02_histogram_equal); title("Image 2 (Histogram equal)"); axis off;
%% Prob 3 : 
image01_ycc = YC_rC_b(image01); image02_ycc = YC_rC_b(image02);
figure(4)
subplot(2,2,1); imagesc(image01); title("Image 1 (Original)"); axis off;
subplot(2,2,2); imagesc(image01_ycc); title("Image 1 (YC_{r}C_{b})"); axis off;
subplot(2,2,3); imagesc(image02); title("Image 2 (Original)"); axis off;
subplot(2,2,4); imagesc(image02_ycc); title("Image 2 (YC_{r}C_{b})"); axis off;
%% Histogram distribution
figure(5)
Image_Histogram_show(image01,image02,"Original image1","Original image2")
figure(6)
Image_Histogram_show(image01_gamma1,image02_gamma1,"image1 (gamma = 0.45)","image2 (gamma = 0.45)")
figure(7)
Image_Histogram_show(image01_gamma2,image02_gamma2,"image1 (gamma = 2.2)","image2 (gamma = 2.2)")
figure(8)
Image_Histogram_show(image01_histogram_equal,image02_histogram_equal,"image1 Histogram Equal","image2 Histogram Equal")
figure(9)
Image_Histogram_show(image01_ycc,image02_ycc,"image1 YC_{r}C_{b}","image2 YC_{r}C_{b}")
%% Functions 
function output_image = Gamma_correction(input_image,gamma)
    image = double(input_image);
    image = image.^(1/gamma);
    image = (image-min(image(:)))/(max(image(:))-min(image(:))) * 255;
    output_image = uint8(image);
end


function Output_image = Histogram_Equalization(input_image)
    R = double(input_image(:,:,1)); G = double(input_image(:,:,2)); B = double(input_image(:,:,3));
    R_distribution = hist(R(:),0:1:255);
    G_distribution = hist(G(:),0:1:255);
    [B_distribution,bin] = hist(B(:),0:1:255);
    Cumculate_R = cumsum(R_distribution);
    Cumculate_G = cumsum(G_distribution);
    Cumculate_B = cumsum(B_distribution);
    Output_image_R = zeros(size(input_image,1),size(input_image,2));
    Output_image_G = zeros(size(input_image,1),size(input_image,2));
    Output_image_B = zeros(size(input_image,1),size(input_image,2));
    for i=1:256
        Output_image_R(R==bin(i)) = (Cumculate_R(i)/Cumculate_R(end)) * 255;
        Output_image_G(G==bin(i)) = (Cumculate_G(i)/Cumculate_G(end)) * 255;
        Output_image_B(B==bin(i)) = (Cumculate_B(i)/Cumculate_B(end)) * 255;
    end
    Output_image_R = (Output_image_R-min(Output_image_R(:)))/(max(Output_image_R(:))-min(Output_image_R(:)))*255;
    Output_image_G = (Output_image_G-min(Output_image_G(:)))/(max(Output_image_G(:))-min(Output_image_G(:)))*255;
    Output_image_B = (Output_image_B-min(Output_image_B(:)))/(max(Output_image_B(:))-min(Output_image_B(:)))*255;
    Output_image = cat(3,Output_image_R,Output_image_G,Output_image_B);
    Output_image = uint8(Output_image);
end


function output_image = YC_rC_b(input_image)
    input_image = double(input_image);
    Transform_matrix_A = [0.299,0.587,0.144;
                          -0.168736,-0.331264,0.5;
                          0.5,-0.418688,-0.081312];
    Transform_matrix_B = [0;128;128];
    middle_image = zeros(size(input_image,1),size(input_image,2),3);
    for i=1:size(input_image,1)
        for j=1:size(input_image,2)
            RGB = reshape( input_image(i,j,:) ,3,1);
            middle_image(i,j,:) = ( Transform_matrix_A * RGB + Transform_matrix_B )'; 
        end
    end
    
    
    R = round(middle_image(:,:,1));  
    [R_distribution,bin] = hist(R(:),min(R(:)):1:max(R(:)));
    Cumculate_R = cumsum(R_distribution); Output_image_R = zeros(size(input_image,1),size(input_image,2));
    for i=1:length(bin)
        Cumculate_R(i) = (Cumculate_R(i)/Cumculate_R(end)) * bin(end);
        Output_image_R(R==bin(i)) = Cumculate_R(i);
    end
    
    
    
    middle_image(:,:,1) = Output_image_R;
    output_image = zeros(size(input_image,1),size(input_image,2),3);
    for i=1:size(input_image,1)
        for j=1:size(input_image,2)
            RGB = reshape(middle_image(i,j,:) ,3,1);
            output_image(i,j,:) = ( Transform_matrix_A \ (RGB - Transform_matrix_B) )'; 
        end
    end
    
    output_image = (output_image-min(output_image(:)))/(max(output_image(:))-min(output_image(:))) * 255;
    output_image = uint8(output_image);
end

function Image_Histogram_show(Image1,Image2,string1,string2)
    [Image_1R,Image_1G,Image_1B] = image_histogram(Image1,0:1:255);
    [Image_2R,Image_2G,Image_2B] = image_histogram(Image2,0:1:255);
    subplot(2,3,1); plot(0:1:255,Image_1R);  title(string1 + "'s Red   distribution")
    subplot(2,3,2); plot(0:1:255,Image_1G);  title(string1 + "'s Green distribution")
    subplot(2,3,3); plot(0:1:255,Image_1B);  title(string1 + "'s Blue  distribution")
    subplot(2,3,4); plot(0:1:255,Image_2R);  title(string2 + "'s Red   distribution")
    subplot(2,3,5); plot(0:1:255,Image_2G);  title(string2 + "'s Green distribution")
    subplot(2,3,6); plot(0:1:255,Image_2B);  title(string2 + "'s Blue  distribution")
end


function [R_distribution,G_distribution,B_distribution] = image_histogram(input_image,bin)
    R = input_image(:,:,1); G = input_image(:,:,2); B = input_image(:,:,3);
    R_distribution = hist(R(:),bin);
    G_distribution = hist(G(:),bin);
    B_distribution = hist(B(:),bin);
end
