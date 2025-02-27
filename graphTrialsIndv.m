%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/21/2025 (MATLAB R2024a)

%this program plots COP and velocity for each trial

%-------------------------------------------------------------------------%

function graphTrialsIndv (trialnames, data, trialInsole,path_Data)

%controls gradient color
colorWay = 'parula'; %choices: 'parula' 'turbo' 'hsv' 'hot' 'cool' 'spring' 'summer' 'autum' 'winter' 'jet'

    for i=1:8
    
        temp = data{i};
    
        leftInsole = temp{1,3}; rightInsole = temp{1,4};
    
        figure(i)
    
        subplot(2,2,1);
        I = imread("INSOLEL.png");
        %imshow(I)
        I = imrotate(I,180); %idk why the left one specifically is weird
        desiredRows = trialInsole.forefootWidth;
        desiredColumns = trialInsole.length;
        K = imresize(I, [desiredColumns, desiredRows]); %scales image to right real world size
        %[rows, columns, numberOfColorChannels] = size(K);% double check for ts
        imshow(K)
        
        hold on
        patch(leftInsole.COP_X, leftInsole.COP_Y,leftInsole.time,'EdgeColor', 'interp');
        colormap(colorWay); %gradient helps match up data to time
        axis on;
        title("Left Insole COP");
        xlabel ("distance (mm)"); ylabel ("distance (mm)");
        hold off
    
        Lavg_MedLat (i)= mean(leftInsole.COP_X, "omitmissing");
        Lavg_AntPost (i)= mean(leftInsole.COP_Y, "omitmissing");
        LmaxDisp_MedLat (i)= max(leftInsole.COP_X)-min(leftInsole.COP_X);
        LmaxDisp_AntPost (i)= max(leftInsole.COP_Y)-min(leftInsole.COP_Y);
    
        subplot(2,2,3);
        I = imread("INSOLER.png");
        desiredRows = trialInsole.forefootWidth;
        desiredColumns = trialInsole.length;
        J = imresize(I, [desiredColumns, desiredRows]);
        %[rows, columns, numberOfColorChannels] = size(K); % double check for ts
        imshow(J)
        
        hold on
        patch(rightInsole.COP_X,rightInsole.COP_Y,rightInsole.time,'EdgeColor', 'interp');
        colormap(colorWay);
        axis on
        title("Right Insole COP");
        xlabel ("distance (mm)"); ylabel ("distance (mm)");
        hold off
    
        Ravg_MedLat (i)= mean(rightInsole.COP_X, "omitmissing");
        Ravg_AntPost (i)= mean(rightInsole.COP_Y, "omitmissing");
        RmaxDisp_MedLat (i)= max(rightInsole.COP_X)-min(rightInsole.COP_X);
        RmaxDisp_AntPost (i)= max(rightInsole.COP_Y)-min(rightInsole.COP_Y);
    
        sgtitle(trialnames(i));
    
    
        %temporal parameters: time to peak COP Velocity
    
        %calculations
        vel_LCOPX = (diff(leftInsole.COP_X)./1000)./diff(leftInsole.time); %mm to m
        vel_LCOPY = (diff(leftInsole.COP_Y)./1000)./diff(leftInsole.time); %mm to m
        vel_LCOP = sqrt(vel_LCOPX.^2 + vel_LCOPY.^2);
        [V_peak_L, peak_idx_L] = max(vel_LCOP);
        L_velPeak (i) = V_peak_L;
        t_peak_L (i)= leftInsole.time(peak_idx_L + 1);
        
    
    
        vel_RCOPX = (diff(rightInsole.COP_X)./1000)./diff(rightInsole.time); %mm to m
        vel_RCOPY = (diff(rightInsole.COP_Y)./1000)./diff(rightInsole.time); %mm to m
        vel_RCOP = sqrt(vel_RCOPX.^2 + vel_RCOPY.^2);
        [V_peak_R, peak_idx_R] = max(vel_RCOP);
        R_velPeak (i) = V_peak_R;
        t_peak_R (i)= rightInsole.time(peak_idx_R + 1);
    
        %scale axis
        ylim_Max = 0;
    
        if(V_peak_L > V_peak_R)
            ylim_Max = V_peak_L;
        else
            ylim_Max = V_peak_R;
        end %end if
    
        %graphs
        subplot(2,2,2);
        patch(leftInsole.time(1:end-1), vel_LCOP, leftInsole.time(1:end-1),'EdgeColor', 'interp'); %time array is trimed b/c od diff function to calc velocity
        colormap(colorWay);
        xlabel("time (s)")
        ylim([0 ylim_Max]); ylabel("velocity m/s");
        title("Left Insole: COP Velocity");
    
        subplot(2,2,4);
        patch(rightInsole.time(1:end-1), vel_RCOP, rightInsole.time(1:end-1),'EdgeColor', 'interp');
        colormap(colorWay);
        xlabel("time (s)");
        ylim([0 ylim_Max]); ylabel("velocity m/s");
        title("Right Insole: COP Velocity");
    
    
        %velocity parameters: mean COP velocity
    
        meanVel_COPL (i) = mean(vel_LCOP,"omitmissing");
        meanVel_COPR (i) = mean(vel_RCOP,"omitmissing");
    
        data_InsoleSpatial = vertcat(Lavg_AntPost,Ravg_AntPost,Lavg_MedLat,Ravg_MedLat, ...
            LmaxDisp_AntPost,RmaxDisp_AntPost,LmaxDisp_MedLat,RmaxDisp_MedLat,...
            meanVel_COPL, meanVel_COPR, t_peak_L,L_velPeak,t_peak_R,R_velPeak);
    
        %save graph 
        figName_= trialnames(i)';
        saveas(figure(i),fullfile(path_Data,figName_), 'jpeg');
    
    end%end for loop

end%end fuunction