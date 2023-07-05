%% Clear everything
clear all; close all; clc;
%% Download image
Giraffe = struct; Grain = struct;
Giraffe.Original = imread('Original_image/giraffe.png');
Grain.Original = imread('Original_image/grain.jpg');
%% Image Segmentation : Watershed Algorithm
Prob1_gap_Threshold = 1;  Prob1_Flooding_level = round(256/Prob1_gap_Threshold)+1;
Giraffe.Watershed = Watershed(Giraffe.Original,Prob1_gap_Threshold,Prob1_Flooding_level);
Giraffe.Watershed_inverse = Watershed(Giraffe.Original,Prob1_gap_Threshold,Prob1_Flooding_level,1);
Grain.Watershed = Watershed(Grain.Original,Prob1_gap_Threshold,Prob1_Flooding_level);
Grain.Watershed_inverse = Watershed(Grain.Original,Prob1_gap_Threshold,Prob1_Flooding_level,1);


figure(1);
subplot(2,3,1);
imagesc(Giraffe.Original); title("Giraffe original")
subplot(2,3,2); colormap gray;
imagesc(Giraffe.Watershed); title("Giraffe watershed")
subplot(2,3,3); colormap gray;
imagesc(Giraffe.Watershed_inverse); title("Giraffe watershed(inverse)")
subplot(2,3,4);
imagesc(Grain.Original); title("Grain original")
subplot(2,3,5); colormap gray;
imagesc(Grain.Watershed); title("Grain watershed")
subplot(2,3,6); colormap gray;
imagesc(Grain.Watershed_inverse); title("Grain watershed(inverse)")


%% Edge Detection : Canny Edge Detection
Prob2_Low_threshold = 0.1; Prob2_High_threshold = 0.4;
Giraffe.Canny = Canny(Giraffe.Original,Prob2_Low_threshold,Prob2_High_threshold);
Prob2_Low_threshold = 0.2; Prob2_High_threshold = 0.63;
Grain.Canny = Canny(Grain.Original,Prob2_Low_threshold,Prob2_High_threshold);

figure(2);
subplot(2,2,1);
imagesc(Giraffe.Original); title("Giraffe original")
subplot(2,2,2); colormap gray;
imagesc(Giraffe.Canny); title("Giraffe Canny")
subplot(2,2,3)
imagesc(Grain.Original); title("Grain original")
subplot(2,2,4); colormap gray;
imagesc(Grain.Canny); title("Grain Canny")

%% Functions for Segmentation : Watershed
function output_image = Watershed(input_image,Threshold,Flooding_level,option)
        if(exist('option')~=0)
            input_image = 255-input_image;
        end
        if(size(input_image,3)==3)
            input_image = rgb2gray(input_image);
        end
        
        Select_value = min(input_image(:));
        output_image = zeros(size(input_image));
        Color = 1;
        for k=1:Flooding_level
            for i=1:size(input_image,1)
                for j=1:size(input_image,2)
                    if(input_image(i,j)>=Select_value && input_image(i,j)<Select_value+Threshold)                    
                        output_image(i,j) = maximum(output_image,i,j);
                        if(output_image(i,j)==0)
                            output_image(i,j) = Color;
                            Color = Color + 1;
                        end
                    end
                end
            end
            Select_value = Select_value+Threshold;
        end
        if(exist('option')~=0)
            output_image = max(output_image(:))-output_image+1;
        end
        
end

function output = maximum(input,x,y)
    if(x==1)
        start_x = x;
        final_x = x+1;
    elseif(x==size(input,1))
        start_x = x-1;
        final_x = x;
    else
        start_x = x-1;
        final_x = x+1;
    end

    if(y==1)
        start_y = y;
        final_y = y+1;
    elseif(y==size(input,2))
        start_y = y-1;
        final_y = y;
    else
        start_y = y-1;
        final_y = y+1;
    end
    output = max(input(start_x:final_x,start_y:final_y),[],"all");

end

%% Functions for Edge detection : Canny
function output_image = Canny(input_image,thresholdRatio_low,thresholdRatio_high)
    if(size(input_image,3)==3)
        input_image = rgb2gray(input_image);
    end
    filter = Gaussian_filter(1,0.4);
    input_image = convn(input_image,filter,"same");
    kx = [-1,0,1;
          -2,0,2;
          -1,0,1];
    ky = [-1,-2,-1;
         0,0,0;
         1,2,1];

     Gx = convn(input_image,kx,"same");
     Gy = convn(input_image,ky,"same");
     
     G = sqrt(Gx.^2+Gy.^2);
     
     theta = rad2deg(atan(Gy./Gx));
     Surpression = non_maximum(G,theta);
     Thresholded_image = signal_threshold(Surpression,thresholdRatio_high,thresholdRatio_low);
     output_image = Edge_tracking(Thresholded_image,25,255);

end 

function filter = Gaussian_filter(filter_size,variance)
    if(exist('sigma')==0)
        variance = 1;
    end
    [x,y] = meshgrid(-filter_size:filter_size,-filter_size:filter_size);
    normalize = 2*pi*variance;
    filter =  exp(-((x.^2 + y.^2) / (2*variance))) / normalize;
end


function output_image = non_maximum(input_image,theta)
    output_image = zeros(size(input_image));
     for i=2:size(input_image,1)-1
        for j=2:size(input_image,2)-1
            if( theta(i,j)>-22.5 && theta(i,j)<22.5 )
                forward = input_image(i,j+1);
                backward= input_image(i,j-1); 
            elseif(theta(i,j)>=22.5 && theta(i,j)<67.5)
                forward = input_image(i+1,j+1);
                backward= input_image(i-1,j-1); 
                
            elseif(theta(i,j)>-67.5 && theta(i,j)<=-22.5)
                forward = input_image(i+1,j-1);
                backward= input_image(i-1,j+1); 
            else
                forward = input_image(i+1,j);
                backward= input_image(i-1,j); 
            end
            
            if( input_image(i,j)>=forward && input_image(i,j)>=backward  )
                output_image(i,j) = input_image(i,j);
            end
        end
     end
end

function output_image = signal_threshold(input_image,thresholdRatio_high,thresholdRatio_low)
     output_image = zeros(size(input_image));
     HighThreshold = thresholdRatio_high*max(input_image(:));
     LowThreshold = thresholdRatio_low*max(input_image(:));
     
     weak_signal = 25;
     strong_signal = 255;
     output_image(input_image>=HighThreshold) = strong_signal;
     output_image(input_image>LowThreshold & input_image<HighThreshold) = weak_signal;
end
    
function output_image = Edge_tracking(input_image,weak_signal,strong_signal)
    output_image = input_image;
    for i=2:size(input_image,1)-1
         for j=2:size(input_image,2)-1
             if(output_image(i,j)==weak_signal)
                if(sum(output_image(i-1:i+1,j-1:j+1)>=strong_signal-1,"all")>0)
                    output_image(i,j) = strong_signal;             
                else
                    output_image(i,j) = 0;
                end
             end
         end
    end
    output_image = output_image>0;
end

     

