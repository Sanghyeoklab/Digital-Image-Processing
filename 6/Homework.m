%% Clear everything
clear all; close all; clc;
%% Load total image
Original_images = struct;
Original_images.Original = double(imread("blur image/hara_original.png"));
Original_images.Cylind_gauss = double(imread("blur image/hara_cylind_gauss.png"));
Original_images.Cylind_blur = double(imread("blur image/hara_cylind_blur.png"));
Original_images.Real_blur = double(imread("blur image/hara_real_blur.png"));
Original_images.Real_gauss = double(imread("blur image/hara_real_gauss.png"));

%% Prob 1
i_T = 10;
filter1_Inverse = Select_filter(Original_images,1,i_T);
plot_images(filter1_Inverse,"Inverse filtered",1)

%% Prob 2
i_T = 10; ratio = 0.99;
filter2_Wiener = Select_filter(Original_images,2,ratio,i_T);
plot_images(filter2_Wiener,"Wiener filtered",2)

%% Prob 3
i_T = 10; gamma = 0.02;
filter3_Constrained_matrix_inversion = Select_filter(Original_images,3,gamma,i_T);
plot_images(filter3_Constrained_matrix_inversion,"CMI filtered",3)

%% Discussion
i_T = 10; gamma = 0.02;
filter4_Constrained_matrix_inversion_v2 = Select_filter(Original_images,4,gamma,i_T);
plot_images(filter4_Constrained_matrix_inversion_v2,"CMI filtered",4)

%% Function
function output_image = inverse_filtered_image(input_image,i_T)
        G = fftn(input_image);
        G1 = real(G); G2 = imag(G);
        [N1,N2] = size(input_image);
        
        pimN = repmat((0:N2-1),N1,1)*pi/N2;
       
        F1 = i_T*sin(pimN).*(  G1.*cos(pimN*(i_T-1)) - G2.*sin(pimN*(i_T-1)) )./sin(pimN*i_T);
        F2 = i_T*sin(pimN).*(  G1.*sin(pimN*(i_T-1)) + G2.*cos(pimN*(i_T-1)) )./sin(pimN*i_T);
        for k=1:i_T-1
            if( mod((k-1)*N2,i_T)==0 )
                m = (k-1)*N2/i_T+1;
                F1(:,m) = G1(:,m);  
                F2(:,m) = G2(:,m);
            end
        end
        
        m0 = floor(N2/i_T);
        F1(:,m0:N2-m0) = G1(:,m0:N2-m0);
        F2(:,m0:N2-m0) = G2(:,m0:N2-m0);
        
        F = F1 + F2*1j;
        output_image = ifftn(F);
end

function output_image = Wiener_filtered_image(input_image,ratio,i_T)
    G = fftn(input_image);
    G1 = real(G); G2 = imag(G);
    N = size(input_image,2);
    pimN = repmat((0:size(input_image,2)-1),size(input_image,1),1)*pi/N;
    
    F1 = i_T*sin(pimN).*sin(pimN*i_T).* ...
        (G1.*cos(pimN*(i_T-1)) - G2.*sin(pimN*(i_T-1)))./ ...
        (sin(pimN*i_T).^2+ratio*(i_T*sin(pimN)).^2  );
    
    F2 = i_T*sin(pimN).*sin(pimN*i_T).* ...
        (G1.*sin(pimN*(i_T-1)) + G2.*cos(pimN*(i_T-1)))./ ...
        (sin(pimN*i_T).^2+ratio*(i_T*sin(pimN)).^2  );
    
    F1(:,1) = G1(:,1)/(1+ratio);  F2(:,1) = G2(:,1)/(1+ratio);
    
    F = F1 + F2*1j;
    output_image = ifftn(F);
end



function output_image = Constrained_matrix_inversion_filtered_image(input_image,ratio,i_T)
    G = fftn(input_image);
    G1 = real(G); G2 = imag(G);
    N = size(input_image,2);
    pimN = repmat((0:size(input_image,2)-1),size(input_image,1),1)*pi/N;
    
    L = 27 -10*cos(2*pimN) -10*cos((N-1)*pimN*2) + 2*cos((N-2)*pimN*2);
    L(:,1) = (N-5)^2 +2 + 2*(N-5)*cos(pimN(:,1)*2) + 2*(N-5)*cos(pimN(:,1)*2*(N-1))...
              + 2*cos(pimN(:,1)*2*(N-2));
    
    
    
    F1 = i_T*sin(pimN).*sin(pimN*i_T).* ...
        (G1.*cos(pimN*(i_T-1)) - G2.*sin(pimN*(i_T-1)))./ ...
        (sin(pimN*i_T).^2+ratio*L.*(i_T*sin(pimN)).^2  );
    
    F2 = i_T*sin(pimN).*sin(pimN*i_T).* ...
        (G1.*sin(pimN*(i_T-1)) + G2.*cos(pimN*(i_T-1)))./ ...
        (sin(pimN*i_T).^2+ratio*L.*(i_T*sin(pimN)).^2  );
    F1(:,1) = G1(:,1);  F2(:,1) = G2(:,1);
    F = F1 + F2*1j;
    output_image = ifftn(F);
end



function output_image = Constrained_matrix_inversion_filtered_image_2(input_image,ratio,i_T)
    G = fftn(input_image);
    G1 = real(G); G2 = imag(G);
    N = size(input_image,2);
    [pimN,pinN] = meshgrid( 0:1:size(input_image,2)-1,0:1:size(input_image,1)-1);
    pinN = pinN*pi/N;
    pimN = pimN*pi/N;
    
    L = abs( -4 + exp(pimN*2*1j) +exp(pimN*2*(N-1)*1j) +exp(pinN*2*(N-1)*1j)+exp(pinN*2)*1j).^2;

    F1 = i_T*sin(pimN).*sin(pimN*i_T).* ...
        (G1.*cos(pimN*(i_T-1)) - G2.*sin(pimN*(i_T-1)))./ ...
        (sin(pimN*i_T).^2+ratio*L.*(i_T*sin(pimN)).^2  );
    
    F2 = i_T*sin(pimN).*sin(pimN*i_T).* ...
        (G1.*sin(pimN*(i_T-1)) + G2.*cos(pimN*(i_T-1)))./ ...
        (sin(pimN*i_T).^2+ratio*L.*(i_T*sin(pimN)).^2  );
    
    F1(:,1) = G1(:,1)./(1+ratio*L(:,1));  F2(:,1) = G2(:,1)./(1+ratio*L(:,1));
    
    F = F1 + F2*1j;
    output_image = ifftn(F);
end





function output_structure = Select_filter(image_structure,filter_number,condition1,condition2)
    output_structure = image_structure;
    if(filter_number==1)
        output_structure.output_image_Cg = inverse_filtered_image(output_structure.Cylind_gauss,condition1);
        output_structure.output_image_Cb = inverse_filtered_image(output_structure.Cylind_blur,condition1);
        output_structure.output_image_Rg = inverse_filtered_image(output_structure.Real_gauss,condition1);
        output_structure.output_image_Rb = inverse_filtered_image(output_structure.Real_blur,condition1);
    elseif(filter_number==2)
        output_structure.output_image_Cg = Wiener_filtered_image(output_structure.Cylind_gauss,condition1,condition2);
        output_structure.output_image_Cb = Wiener_filtered_image(output_structure.Cylind_blur,condition1,condition2);
        output_structure.output_image_Rg = Wiener_filtered_image(output_structure.Real_gauss,condition1,condition2);
        output_structure.output_image_Rb = Wiener_filtered_image(output_structure.Real_blur,condition1,condition2);
    elseif(filter_number==3)
        output_structure.output_image_Cg = Constrained_matrix_inversion_filtered_image(output_structure.Cylind_gauss,condition1,condition2);
        output_structure.output_image_Cb = Constrained_matrix_inversion_filtered_image(output_structure.Cylind_blur,condition1,condition2);
        output_structure.output_image_Rg = Constrained_matrix_inversion_filtered_image(output_structure.Real_gauss,condition1,condition2);
        output_structure.output_image_Rb = Constrained_matrix_inversion_filtered_image(output_structure.Real_blur,condition1,condition2);
    elseif(filter_number==4)
        output_structure.output_image_Cg = Constrained_matrix_inversion_filtered_image_2(output_structure.Cylind_gauss,condition1,condition2);
        output_structure.output_image_Cb = Constrained_matrix_inversion_filtered_image_2(output_structure.Cylind_blur,condition1,condition2);
        output_structure.output_image_Rg = Constrained_matrix_inversion_filtered_image_2(output_structure.Real_gauss,condition1,condition2);
        output_structure.output_image_Rb = Constrained_matrix_inversion_filtered_image_2(output_structure.Real_blur,condition1,condition2);
    end
end


function plot_images(image_structure,label,figure_num)
    figure(figure_num); colormap gray;
    subplot(4,3,1); imagesc(abs(image_structure.Original))
    title("Original image");axis off;
    subplot(4,3,2); imagesc(abs(image_structure.Cylind_gauss))
    title("Cylind gauss image");axis off;
    subplot(4,3,3); imagesc(abs(image_structure.output_image_Cg))
    title(label + " Cylind gauss image");axis off;

    subplot(4,3,4); imagesc(abs(image_structure.Original))
    title("Original image");axis off;
    subplot(4,3,5); imagesc(abs(image_structure.Cylind_blur))
    title("Cylind blur image");axis off;
    subplot(4,3,6); imagesc(abs(image_structure.output_image_Cb))
    title(label + " Cylind blur image");axis off;
    
    subplot(4,3,7); imagesc(abs(image_structure.Original))
    title("Original image");axis off;
    subplot(4,3,8); imagesc(abs(image_structure.Real_gauss))
    title("Real gauss image");axis off;
    subplot(4,3,9); imagesc(abs(image_structure.output_image_Rg))
    title(label + " Real gauss image");axis off;
    
    subplot(4,3,10); imagesc(abs(image_structure.Original))
    title("Original image");axis off;
    subplot(4,3,11); imagesc(abs(image_structure.Real_blur))
    title("Real blur image");axis off;
    subplot(4,3,12); imagesc(abs(image_structure.output_image_Rb))
    title(label + " Real blur image");axis off;
end