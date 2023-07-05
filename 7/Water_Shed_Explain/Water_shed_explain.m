%% Clear everything
clear all; close all;clc;

%% Make f
f = [6, 10, 10, 10, 7, 8, 8, 7, 4, 11, 6, 6, 6, 8, 8, 5, 5, 5, 5, 6, 4, 5, 4];

f = round(rand(size(f))*10)+2;

%% Show each flood
figure(1)
for Flood=1:9
    minimum_value = (min(f(:))+Flood)*ones(size(f));
    color = Watershed(f,Flood);
    subplot(3,3,Flood)
    plot_image(f,minimum_value,color); 
    if(Flood==1)
        title("All minima pierced; flood begins")
    else
        title("Flood level : " + string(Flood-1))
    end
end





function [output,label] = for_plot(input)
    label = zeros(1,numel(input)*2-2); output = zeros(1,numel(input)*2);
    label(1:2:end) = 2:(numel(input));
    label(2:2:end) = 2:(numel(input));
    label = [1,label,numel(input)+1];
    output(1:2:end) = input;
    output(2:2:end) = input;
end
function plot_image(input1,input2,color)

    output = for_plot(input1);
    color = for_plot(color); color = color/9;
    [minimum_value,label] = for_plot(input2);
    hold on;grid on;
    plot(label,output,'k-','LineWidth',1.5)
    plot(label,minimum_value,'k--','LineWidth',1.5)
    
    yticks(0:1:16);yticklabels([]);
    xticks(1:numel(input1)+1);xticklabels([]);
    axis([1,numel(input1)+1,0,16]);
    for i=1:numel(color)-1
        if( output(i)<minimum_value(i) && output(i+1)<minimum_value(i+1) )
            x = [label(i),label(i),label(i+1),label(i+1)];
            y = [output(i),minimum_value(i),minimum_value(i+1),output(i+1)];
            fill(x,y,[color(i),color(i),color(i)],'EdgeColor','none');
        end
    end
    
end


function output_image = Watershed(input_image,Flooding_level)
  
        Threshold_value = min(input_image(:));
        output_image = zeros(size(input_image));
        index = 1;
        for k=1:Flooding_level
        
            for i=1:size(input_image,2)
                
                    if(input_image(i)==Threshold_value)
                        output_image(i) = maximum(output_image,i);
                        if(output_image(i)==0)
                            output_image(i) = index;
                            index = index + 1;
                        end
                    end

            end
            Threshold_value = Threshold_value+1;
        end
end

function output = maximum(input,x)
    if(x==1)
        start_x = x;
        final_x = x+1;
    elseif(x==size(input,2))
        start_x = x-1;
        final_x = x;
    else
        start_x = x-1;
        final_x = x+1;
    end
    output = max(input(1,start_x:final_x),[],"all");
end
