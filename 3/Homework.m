%% Data clear
clear all; close all; clc % Clear everything.
%% Prob 1
Original_image = imread('Lenna.png'); 
% Read Original image "Lenna.png".
Original_image = rgb2gray(Original_image); 
% Original image is color, so it need to convert gray image.
% Consider 8 bits grey label image which size is 256 X 256.
%% Prob 2
g = double(Original_image);
% For calculate, make int form values to double form values.
g_star = zeros(size(g)); g_star(1,:) = g(1,:); g_star(:,1) = g(:,1);
% Remain First row and First column.
for j=2:256
    for i=2:255
        g_star(i,j) = 1/4*( g_star(i-1,j) + g_star(i-1,j-1) + g_star(i,j-1) + g_star(i+1,j-1));
        % fill by raster scan order.
    end
    g_star(256,j) = 1/3*( g_star(255,j-1) + g_star(255,j) + g_star(256,j-1) );
    % In last row, consider number of value.
end
e = g - g_star;
% Get error of g_star and g.
%% Prob 3
e_flat = reshape(e,1,numel(e));
% Make Image array to 1D array for Calculate
delta_bin = 0.5;
% Unit of error
[PDF,nbin] = hist(e_flat,-255:delta_bin:255); PDF = PDF/(sum(PDF)*delta_bin);
% Make PDF by histogram. Consider PDF is probability for range so divide by
% unit of error.
figure(1)
plot(nbin,PDF); title("PDF of e"); xlabel("e"); ylabel("Probability")
% Show PDF.
d = prctile(unique(e_flat),0:100/64:100); 
% Set initial d by percentage values of unique values.
d(1) = -255-1; d(end) = 255+1; d_after = d; 
% Consider first and end of d is edge of bin.
while(1)
    for i=1:64
        Range = ( e_flat > d(i) ) & ( e_flat < d(i+1) );
        % Adopt only values ??in each range.
        r(i) = mean(e_flat(Range));
        % Actually it has to use PDF but anyway take mean value of adopted
        % values
    end
    for i=1:63
        d_after(i+1) = (r(i) + r(i+1))/2;
        % Make d average by neighbor r 
    end
    if(isequal(round(d,3),round(d_after,3)))
    % Do this repete until d become convergence.
        signal = zeros(size(e_flat));
        % Match e with correspond symbol. Actually I make function for only
        % natural values.
        for i=1:64
            Range = ( e_flat > d(i) ) & ( e_flat < d(i+1) );
            % Get Probability of each symbol.
            Probability(i,1) = sum(Range);
            signal(Range) = i;
        end
        Probability = Probability/sum(Probability);
        break;
    else
        d = d_after;
    end
end
figure(2)
plot(r,Probability); xlabel("Symbol"); ylabel("Probability");
title("Probability of each symbol")

%% Prob 4
% dict = huffmandict(sort(unique(signal)),Probability); % This is Matlab code
% Matlab given function and function that I made is not compatible!!!!!!!
Dictionary = Huffman_Dict(Probability);
% Make dictionary by function that I made.
fprintf("==================================================================\n")
fprintf("Table\n");
fprintf("%-3s (%7s ~ %-7s)  -> %8s\t%-8s\t\t%s\n","i","d(i)","d(i+1)","r(i)","Probability","Dictionary")
fprintf("==================================================================\n")
for i=1:64
    fprintf("%-3d (%7.2f ~ %-7.2f)  -> %9.2f\t%8.4f\t\t%s\n",i,d(i),d(i+1),r(i),Probability(i),Dictionary(i))
end
fprintf("==================================================================\n\n\n\n")

x(1) = d(1); x(2:2:127) = d(2:end-1); x(3:2:127) = d(2:end-1); x(128) = d(end);
y(1:2:128) = r; y(2:2:128) = r;
figure(3)
plot(x,y); xlabel("e"); ylabel("e^{*}"); title("Typical input-output characteristic of a quantizer")
grid on;

%% Prob 5
% hcode = huffmanenco(signal,dict); % This is Matlab code
% Matlab given function and function that I made is not compatible!!!!!!!
hcode = Huffman_encode(signal,Dictionary);
% Make code by function that I made.
entropy = sum(-Probability.*log2(Probability));
fprintf("Entropy : %f\n",entropy);
fprintf("Length of Code : %d => %f per pixel\n",length(hcode),length(hcode)/numel(e_flat))
% Show length of code that I encoding by huffman code.
%% Prob 6
% dsignal = huffmandeco(hcode,dict); % This is Matlab code
% Matlab given function and function that I made is not compatible!!!!!!!
dsignal = Huffman_decode(hcode,Dictionary);

e_star = zeros(size(dsignal));
for i=1:64
    e_star(dsignal==i) = r(i);
end
% Reconstruct error by correspond quantization value.
e_star = reshape(e_star,256,256);
g_tilde = zeros(size(g));
g_tilde(1,:) = g(1,:); g_tilde(:,1) = g(:,1);
% Consider remain first row and column. So their value never change.
for j=2:256
    for i=2:255
        g_tilde(i,j) = 1/4*( g_star(i-1,j) + g_star(i-1,j-1) + g_star(i,j-1) + g_star(i+1,j-1))...
                      + e_star(i,j);
    end
    g_tilde(256,j) = 1/3*( g_star(255,j-1) + g_star(255,j) + g_star(256,j-1) )...
                    + e_star(i,j);
end
% Reconstruct by DPCM.
%% Prob 7

% Show Original image and reconstruct image. Also show error map with PSNR.
MAXI = max(g(:));
MSE = mean( (g(:) - g_tilde(:)).^2 );
PSNR = 20*log10(MAXI) - 10*log10(MSE);
fprintf("PSNR : %f dB\n",PSNR);
figure(4)
colormap gray;
subplot(2,2,1); imagesc(Original_image); title("\fontsize{16}Original image");
subplot(2,2,2); imagesc(g_tilde);  title('\fontsize{16},$\tilde{g}$  (Code length : ' + string(length(hcode)) + ')','Interpreter','latex')

subplot(2,2,[3,4]); imagesc(g - g_tilde); 
title("\fontsize{16}Error map    PSNR : " + string(PSNR));
colorbar;
%% Functions for Huffman Code.
% Function for make Huffman dictionary
function Dictionary = Huffman_Dict(Probability)
    Dic_size = length(Probability);
    % Initialize Dictionary.
    for i=1:Dic_size
        Dictionary(i,1) = "";
    end
    
    % For contain index.
    index = zeros(Dic_size,Dic_size);
    index(:,1) = 1:Dic_size;

    
    for i=1:Dic_size-1
        [Probability,Sort_index] = sort(Probability,'descend');
        index = index(Sort_index,:);
        % Each iteration, sort by descend.
        for j = index(end,:)
            if(j==0)
                break
            end
            Dictionary(j) = "1" + Dictionary(j) ;
        end

        for j = index(end-1,:)
            if(j==0)
                break
            end
            Dictionary(j) = "0" + Dictionary(j);
        end
        % Last row that have low probability, add signal in each condition.
        Probability(end-1) = Probability(end-1) + Probability(end); Probability(end) = [];
        % Combine last two row.
        start = find(index(end-1,:)==0,1); index_num = find(index(end,:)~=0);
        index(end-1,start:start+length(index_num)-1) = index(end,index_num);
        index(end,:) = [];
    end
end

% Function for make Huffman code, symbol is natural number.
function Code = Huffman_encode(symbol,Dictionary)
    Code = "";
    for i=1:length(symbol)
        Code = Code + Dictionary(symbol(i));
    end
    Code = char(Code);
end

% Function for decode Huffman code.
function Recon = Huffman_decode(Code,Dictionary)
    Recon = [];
    while(length(Code)~=0)
        for i=1:length(Dictionary)
            if(length(Code)>=length( char(Dictionary(i)) ))
            if(Code(1:length( char(Dictionary(i)) )) == Dictionary(i))
                Recon(end+1) = i;
                Code(1:length(char(Dictionary(i)))) = [];
                break;
            end
            end
        end
    end
    Recon = char(Recon);
end