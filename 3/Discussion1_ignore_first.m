%% Data clear
% If you want to see more detail, watch Homework code. It hust for
% discussion that certain case.
clear all; close all; clc;
%% Prob 1
Original_image = imread('Lenna.png'); 
Original_image = rgb2gray(Original_image); 
%% Prob 2
g = double(Original_image);
g_star = zeros(size(g)); g_star(1,:) = g(1,:); g_star(:,1) = g(:,1);
for j=2:256
    for i=2:255
        g_star(i,j) = 1/4*( g_star(i-1,j) + g_star(i-1,j-1) + g_star(i,j-1) + g_star(i+1,j-1));
    end
    g_star(256,j) = 1/3*( g_star(255,j-1) + g_star(255,j) + g_star(256,j-1) );
end
% First row and column is same anyway, so ignore them in calculate.
e = g(2:end,2:end) - g_star(2:end,2:end);
%% Prob 3
e_flat = reshape(e,1,numel(e));
delta_bin = 0.5;
[PDF,nbin] = hist(e_flat,-255:delta_bin:255); PDF = PDF/(sum(PDF)*delta_bin);
figure(1)
plot(nbin,PDF); title("PDF of e"); xlabel("e"); ylabel("Probability")
d = prctile(unique(e_flat),0:100/64:100); 
d(1) = -255-1; d(end) = 255+1; d_after = d; 
while(1)
    for i=1:64
        Range = ( e_flat > d(i) ) & ( e_flat < d(i+1) );
        r(i) = mean(e_flat(Range));
    end
    for i=1:63
        d_after(i+1) = (r(i) + r(i+1))/2;
    end
    if(isequal(round(d,3),round(d_after,3)))
        signal = zeros(size(e_flat));
        for i=1:64
            Range = ( e_flat > d(i) ) & ( e_flat < d(i+1) );
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
Dictionary = Huffman_Dict(Probability);
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
hcode = Huffman_encode(signal,Dictionary);
entropy = sum(-Probability.*log2(Probability));
fprintf("Entropy : %f\n",entropy);
fprintf("Length of Code : %d => %f per pixel\n",length(hcode),length(hcode)/numel(e_flat))
%% Prob 6
dsignal = Huffman_decode(hcode,Dictionary);
e_star = zeros(size(dsignal));
for i=1:64
    e_star(dsignal==i) = r(i);
end

e_star = reshape(e_star,255,255);
g_tilde = zeros(size(g));
g_tilde(1,:) = g(1,:); g_tilde(:,1) = g(:,1);
for j=2:256
    for i=2:255
        g_tilde(i,j) = 1/4*( g_star(i-1,j) + g_star(i-1,j-1) + g_star(i,j-1) + g_star(i+1,j-1))...
                      + e_star(i-1,j-1);
    end
    g_tilde(256,j) = 1/3*( g_star(255,j-1) + g_star(255,j) + g_star(256,j-1) )...
                    + e_star(i-1,j-1);
end
%% Prob 7

MAXI = max(g(2:end));
MSE = mean( (g(2:end) - g_tilde(2:end)).^2 );
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
function Dictionary = Huffman_Dict(Probability)
    Dic_size = length(Probability);
    for i=1:Dic_size
        Dictionary(i,1) = "";
    end
    
    index = zeros(Dic_size,Dic_size);
    index(:,1) = 1:Dic_size;

    
    for i=1:Dic_size-1
        [Probability,Sort_index] = sort(Probability,'descend');
        index = index(Sort_index,:);
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
        Probability(end-1) = Probability(end-1) + Probability(end); Probability(end) = [];
        start = find(index(end-1,:)==0,1); index_num = find(index(end,:)~=0);
        index(end-1,start:start+length(index_num)-1) = index(end,index_num);
        index(end,:) = [];
    end
end

function Code = Huffman_encode(symbol,Dictionary)
    Code = "";
    for i=1:length(symbol)
        Code = Code + Dictionary(symbol(i));
    end
    Code = char(Code);
end

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