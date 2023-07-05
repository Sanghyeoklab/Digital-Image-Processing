%% Clear everything
clear all; close all; clc;
folder = "Image_files\";
Orginal_image = double(imread(folder+"crowd_original.bmp"));
image_gaussian_noise = double(imread(folder+"crowd_gau_30%.bmp"));
image_SP_noise = double(imread(folder+"crowd_sp_30%.bmp"));
% load each images.

%% Prob 1
Gaussian_median = median_filtered_image(image_gaussian_noise);
image_SP_median = median_filtered_image(image_SP_noise);
% Apply median filter in each images.
figure(1); colormap gray;
subplot(2,3,1); imagesc(image_gaussian_noise); title("30% Gaussian noise")
subplot(2,3,2); imagesc(Gaussian_median); title("Apply meidan filter in Gaussian noise")
subplot(2,3,3); imagesc(Orginal_image - Gaussian_median); title("Gaussian Error map   MAE : " + string(mean(abs(Orginal_image(:) - Gaussian_median(:)))))
subplot(2,3,4); imagesc(image_SP_noise); title("30% Salt & Pepper noise")
subplot(2,3,5); imagesc(image_SP_median); title("Apply meidan filter in Salt & Pepper noise")
subplot(2,3,6); imagesc(Orginal_image - image_SP_median); title("Salt & Pepper Error map   MAE : " + string(mean(abs(Orginal_image(:) - image_SP_median(:)))))
% Show original error image, filtered image, error map with original image.
% & MAE 
%% Prob 2
Hamming_filter_5_5 = Hamming_filter(5,pi/2);
% Get 5 x 5 filter by hamming window.
Gaussian_Hamming = Hamming_filtered_image(image_gaussian_noise,Hamming_filter_5_5);
image_SP_Hamming = Hamming_filtered_image(image_SP_noise,Hamming_filter_5_5);
% Apply filter.
figure(2); colormap gray;
subplot(2,3,1); imagesc(image_gaussian_noise); title("30% Gaussian noise")
subplot(2,3,2); imagesc(Gaussian_Hamming); title("Apply Hamming filter in Gaussian noise")
subplot(2,3,3); imagesc(Orginal_image - Gaussian_Hamming); title("Gaussian Error map   MAE : " + string(mean(abs(Orginal_image(:) - Gaussian_Hamming(:)))))
subplot(2,3,4); imagesc(image_SP_noise); title("30% Salt & Pepper noise")
subplot(2,3,5); imagesc(image_SP_Hamming); title("Apply Hamming filter in Salt & Pepper noise")
subplot(2,3,6); imagesc(Orginal_image - image_SP_Hamming); title("Salt & Pepper Error map   MAE : " + string(mean(abs(Orginal_image(:) - image_SP_Hamming(:)))))
%% Discussion 1 : Gaussian noise 
Gaussian_error = Orginal_image - double(image_gaussian_noise);
% Get gaussian error by subtraction.
Orginal_image_fft = fftshift((fftn(Orginal_image)));
% Get Orignal image's frequency domain.
Gaussian_error_fft = fftshift((fftn(Gaussian_error)));
% Get frequency domain of gaussian error by subtraction.
Gaussian_error_fft_remove = Gaussian_error_fft;
Gaussian_error_fft_remove(256-2:256+2,256-2:256+2) = 0;
% For remove center region for watch other region.
[Gaussian_error_histogram,nbin] = hist(Gaussian_error(:),-256*0.3:1:256*0.3);
Gaussian_error_histogram = Gaussian_error_histogram/sum(Gaussian_error_histogram(:));
% For get error's distribution.
figure(3)
subplot(2,2,1); imagesc(abs(Gaussian_error)); title("Error between Original and 30% Gaussian")
subplot(2,2,2); imagesc(abs(Orginal_image_fft)); title("Fourier Transform of Original image")
subplot(2,2,3); imagesc(abs(Gaussian_error_fft)); title("Fourier Transform of Error")
subplot(2,2,4); imagesc(abs(Gaussian_error_fft_remove)); title("Fourier Transform of Error (For show detail, remove center small patch)")
colormap gray;

figure(4)
plot(nbin,Gaussian_error_histogram); title("Gaussian Error distribution");
xlabel("Error"); ylabel("PMF");
% Plot in each correspond images.
%% Discussion 2 : Hamming Filter
fprintf("Hamming filter : \n")
Hamming_filter_5_5
% For show 5 x 5 filter.
Filter_5_5_extend = extend_signal(Hamming_filter_5_5,513);
% Other signal = 0, just for watch more detail in frequency domain.
Hamming_filter_21_21 = Hamming_filter(21,pi/2);
% Watch more detail about filter with hamming window.
Filter_21_21_extend = extend_signal(Hamming_filter_21_21,513);

figure(5)
subplot(1,2,1)
[w1,w2] = meshgrid(0:1:size(Hamming_filter_5_5,1)-1,0:1:size(Hamming_filter_5_5,2)-1);
mesh(w1,w2,Hamming_filter_5_5); title("5 x 5 Hamming filter")
xlabel("n_{1}"); ylabel("n_{2}");
subplot(1,2,2)
[w1,w2] = meshgrid(0:1:size(Hamming_filter_21_21,1)-1,0:1:size(Hamming_filter_21_21,2)-1);
mesh(w1,w2,Hamming_filter_21_21); title("21 x 21 Hamming filter")
xlabel("n_{1}"); ylabel("n_{2}");
% Plot 5 x 5 & 21 x 21 filter.

figure(6)
subplot(3,2,1)
signalfrequency_show(Filter_5_5_extend); title("FFT of 5 x 5 Hamming Filter")
subplot(3,2,2)
signalfrequency_show_dB_scale(Filter_5_5_extend);title("FFT of 5 x 5 Hamming Filter (dB scale)")
subplot(3,2,3)
signalfrequency_show(Filter_21_21_extend); title("FFT of 21 x 21 Hamming Filter")
subplot(3,2,4)
signalfrequency_show_dB_scale(Filter_21_21_extend);title("FFT of 21 x 21 Hamming Filter (dB scale)")
subplot(3,2,5)
signalfrequency_show_dB_scale(Filter_21_21_extend);title("FFT of 21 x 21 Hamming Filter (dB scale, w1 side watching)")
subplot(3,2,6)
signalfrequency_show_dB_scale(Filter_21_21_extend);title("FFT of 21 x 21 Hamming Filter (dB scale, w2 side watching)")
% Watch in frequency domain.
%% Discusssion 3 : Gaussian filter
Gaussian_filtered_image = Gaussian_filter(image_gaussian_noise,4);
% Apply gaussin filter
figure(7)
colormap gray;
subplot(2,2,1); imagesc(Gaussian_Hamming); title("Apply Hamming filter in Gaussian noise")
subplot(2,2,2); imagesc(Orginal_image - Gaussian_Hamming); title("Hamming filtered image Error map   MAE : " + string(mean(abs(Orginal_image(:) - Gaussian_Hamming(:)))))
subplot(2,2,3); imagesc(Gaussian_filtered_image); title("Apply Gaussian filter in Gaussian noise")
subplot(2,2,4); imagesc(Orginal_image - Gaussian_filtered_image); title("Gaussian filtered image Error map   MAE : " + string(mean(abs(Orginal_image(:) - Gaussian_filtered_image(:)))))
%% Functions for make filter

% For median filter
function Output_image = median_filtered_image(input_image)
    Padding_image = zeros(size(input_image,1)+4,size(input_image,2)+4);
    Padding_image(3:end-2,3:end-2) = input_image;
    % For fit size of input image and output image, two zero - padding in
    % each direction.
    Output_image = zeros(size(input_image));
    for i=1:size(input_image,1)
        for j=1:size(input_image,2)
            Output_image(i,j) = middle_point(Padding_image(i:i+4,j:j+4));
            % Apply middle rank in correspond values.
        end
    end
end

% Extend signal for watch detail in frequency domain. 
function output_signal = extend_signal(input_signal,N)
    output_signal = zeros(N,N);
    output_signal(1:size(input_signal,1),1:size(input_signal,2)) = input_signal;
end

% Plot 3D image.
function signalfrequency_show(signal)
    [w1,w2] = meshgrid(linspace(-pi,pi,size(signal,1)),linspace(-pi,pi,size(signal,2)));
    Frequency = fftshift(fftn(signal));
    mesh(w1,w2,abs(Frequency))
    xticks([-pi -pi/2 0 pi/2 pi]); xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    yticks([-pi -pi/2 0 pi/2 pi]); yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    xlabel("w_{1}"); ylabel("w_{2}")
    colorbar;
end

% Watch 3D image that dB scale.
function signalfrequency_show_dB_scale(signal)
    [w1,w2] = meshgrid(linspace(-pi,pi,size(signal,1)),linspace(-pi,pi,size(signal,2)));
    Frequency = fftshift(fftn(signal));
    mesh(w1,w2,10*log10(abs(Frequency)))
    xticks([-pi -pi/2 0 pi/2 pi]); xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    yticks([-pi -pi/2 0 pi/2 pi]); yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    xlabel("w_{1}"); ylabel("w_{2}")
    colorbar;
end

% Make hamming windowed filter.
function Filter = Hamming_filter(filter_size,cut_off)
    alpha = 0.46; M = filter_size;
    Hamming_window = (1-alpha) - alpha*cos((2*pi*(0:1:M-1))/(M-1));
    % Make hamming window.
    Filter = Hamming_window.*mysinc((0:1:M-1)-(M-1)/2,pi/cut_off);
    % Apply hamming window in correspond sinc function.
    Filter = Filter'*Filter;
    Filter = Filter/sum(Filter(:));
    % Make sum of filter weight 1.
end

function Output_image = Hamming_filtered_image(input_image,Filter)
    Filter_size = size(Filter,1);
    % Get Hamming window
    Padding_image = zeros(size(input_image,1)+Filter_size-1, size(input_image,2)+Filter_size-1);
    Padding_image(1+(Filter_size-1)/2:end-(Filter_size-1)/2,1+(Filter_size-1)/2:end-(Filter_size-1)/2) = input_image;
    % Apply zero padding
    Output_image = zeros(size(input_image));
    for i=1:size(input_image,1)
        for j=1:size(input_image,2)
            Output_image(i,j) = sum( Padding_image(i:i+Filter_size-1,j:j+Filter_size-1).*Filter,"all");
            % Apply filter.
        end
    end
end

% Function for gaussian filter
function [Output_image, filter] = Gaussian_filter(input_image,sigma)
    filter = exp(-(-2:1:2).^2/(2*sigma)); filter = filter'*filter;
    filter = filter/sum(filter(:));
    
    Padding_image = zeros(size(input_image,1)+4, size(input_image,2)+4);
    Padding_image(3:end-2,3:end-2) = input_image;
    Output_image = zeros(size(input_image));
    
    for i=1:size(input_image,1)
        for j=1:size(input_image,2)
            Output_image(i,j) = sum( Padding_image(i:i+4,j:j+4).*filter,"all");
        end
    end
end

% Function for correspond sinc function.
function y = mysinc(x,N)
     y = 1/N * sin(pi*x/N)./((pi*x/N));
     y(isnan(y)) = 1 / N;
end

% Function for get middle rank value.
function point = middle_point(data)
    data = reshape(data,1,numel(data));
    data = sort(data);
    point = data(round(numel(data)/2));
end