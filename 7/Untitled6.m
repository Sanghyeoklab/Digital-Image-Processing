%% Clear everything
clear all; close all; clc;
%% Download image
Giraffe = imread('giraffe.png');
Grain = imread('grain.jpg');
%% Histogramming and Thresholding
input_image = double( rgb2gray(Giraffe) );
% [Average,Variance] = GMM(input_image,0.1,5);
bin = 0:1:255;
distriubution = hist(input_image(:),bin);
distriubution = distriubution/sum(distriubution);

figure(1)
hold on;
plot(bin,distriubution)

% for i=1:length(Average)
%     plot(bin,Gaussian(bin,Average(i),Variance(i)));
% end

% scatter(bin(local),distriubution(local))
% [mask,t] = histogram_threshold(Giraffe);
% figure(1)
% subplot(4,1,1);imagesc(Giraffe);
% subplot(4,1,2);imagesc(mask(1)); colormap gray;
% subplot(4,1,3);imagesc(mask(2)); colormap gray;
% subplot(4,1,4);imagesc(mask(3)); colormap gray;


function [output,t] = histogram_threshold(input_image)
    if(size(input_image,3)==3)
        input_image = rgb2gray(input_image);
    end
    [Average,Variance] = GMM(input_image,0.01);
    t = Threshold(Average,Variance);
    t = [0;t(:);255];
    output = zeros(0,size(input_image,1),size(input_image,2));
    for i=1:length(t)-1
        mask = zeros(1,size(input_image,1),size(input_image,2));
        mask(1,input_image<=t(i)&input_image>=t(i+1))=1;
        output = cat(1,output,mask);
    end
end

function t = Threshold(Average,Variance)
    syms x real;
    for i=1:length(Average)-1
        first = 1/sqrt(2*pi*Variance(i))*exp(-(x-Average(i)).^2/(Variance(i)*2*pi));
        second = 1/sqrt(2*pi*Variance(i+1))*exp(-(x-Average(i+1)).^2/(Variance(i+1)*2*pi));
    end
    x = solve([first==second],x);
    t = double(x);
end

function logic = is_localmax(graph,index,gap)
    if(index<=gap)
        start = 1;
    else
        start = index - gap;
    end
    if(index>=length(graph)-gap)
        final = length(graph);
    else
        final = index + gap;
    end
    
    if( graph(index)>=max(graph(start:final),[],"all") && max(graph(start:final),[],"all")>0 )
        logic = 1==1;
    else
        logic = 1==0;
    end
end




function output = Gaussian(bin,Average,variance)
    output = 1/sqrt(2*pi*variance)*exp(-(bin-Average).^2/(variance*2));
end

function [Average,Variance] = GMM(input_image,X,gap)
    if(~exist('X'))
        gap = 0.1;
    end
    
    if(~exist('gap'))
        gap = 20;
    end
    bin = 0:1:255;
    histogram_distribution = hist(double(input_image(:)),bin);
    histogram_distribution = histogram_distribution/sum(histogram_distribution);
    
    logic = [];
    for i=1:length(histogram_distribution)
        if(is_localmax(histogram_distribution,i,gap))
            logic = [logic,i];
        end
    end
    
    Average = bin(logic)-1; Variance = ones(size(Average))*1;
    
    
    condition = 1;
    while(condition)
        condition = 0;
        Probability = zeros(length(Average),size(input_image,1),size(input_image,2));
        for i=1:length(Average)
            Gaussian_distribution = Gaussian(bin,Average(i),Variance(i));
            for j=1:256
                Probability(i,input_image==bin(j)) = Gaussian_distribution(j);
            end
        end
        Probability = Probability./sum(Probability,1);
        Probability(isnan(Probability)) = 0;
        for i=1:length(Average)
            Probability_i = squeeze(Probability(1,:,:));
            Average_after = sum(Probability_i(:).*input_image(:))/sum(Probability_i(:));
            Variance_after = sum(Probability_i(:).*(input_image(:)-Average_after).^2)/sum(Probability_i(:));
            if( abs(Average(i)-Average_after)>=abs(Average(i))*X )
                condition = 1;
                Average(i) = Average_after;
            end
            if( abs(Variance(i)-Variance_after)>=abs(Variance(i))*X )
                condition = 1;
                Variance(i) = Variance_after;
            end
        end
    end
    
end





function output = Ostu(input)
    distribution = hist(double(input(:)),0:255);
    distribution = distribution/sum(distribution(:));
    P1 = cumsum(distribution);
    P2 = 1-P1;
    m = cumsum( (0:255).*distribution );
    squared_sigma_B = ((m - m(end)*P1).^2)./(P1.*P2);
    [temp, index] = max(squared_sigma_B);
    output = imbinarize(input, (index-1)/255);
end





