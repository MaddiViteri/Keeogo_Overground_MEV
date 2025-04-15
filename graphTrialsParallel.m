%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/18/2025 (MATLAB R2024a)

%this program plots and compares COP trials with similar conditions

%-------------------------------------------------------------------------%

function graphTrialsParallel (trialnames, data, trialInsole, path_Data)

%controls gradient color
colorWayNON = 'b'; colorWayKON = 'm';

%controls opacity
a = 1;

%SS speed  

figure(1)
    %left insole

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

        j=1;
        temp = data{j}; leftInsole1 = temp{1,3};
        h1 = scatter(leftInsole1.COP_X, leftInsole1.COP_Y,colorWayNON,".");

        j = 2;
        temp = data{j}; leftInsole2 = temp{1,3}; 
        h2 = scatter(leftInsole2.COP_X, leftInsole2.COP_Y,colorWayNON,".");

        leftTotalX = [leftInsole1.COP_X; leftInsole2.COP_X];
        leftTotalY = [leftInsole1.COP_Y; leftInsole2.COP_Y];

        p = polyfit(leftTotalX,leftTotalY,2);
        x_fit = linspace(min(leftTotalX),max(leftTotalX),100);
        y_fit = polyval(p,x_fit);

        plot(x_fit,y_fit, 'r');

        axis on;
        title("Left Insole COP NON");
        xlabel ("distance (mm)"); ylabel ("distance (mm)");
        hold off

        % subplot(2,2,3);
        % I = imread("INSOLEL.png");
        % %imshow(I)
        % I = imrotate(I,180); %idk why the left one specifically is weird
        % desiredRows = trialInsole.forefootWidth;
        % desiredColumns = trialInsole.length;
        % K = imresize(I, [desiredColumns, desiredRows]); %scales image to right real world size
        % %[rows, columns, numberOfColorChannels] = size(K);% double check for ts
        % imshow(K)
        % 
        % hold on 
        % 
        % j = 5;
        % temp = data{j}; leftInsole = temp{1,3}; 
        % 
        % h3 = scatter(leftInsole.COP_X, leftInsole.COP_Y,colorWayKON,"."); alpha(a);
        % 
        % j = 6;
        % temp = data{j}; leftInsole = temp{1,3}; 
        % h4 = scatter(leftInsole.COP_X, leftInsole.COP_Y,colorWayKON,"."); alpha(a);
        % 
        % axis on;
        % title("Left Insole COP KON");
        % xlabel ("distance (mm)"); ylabel ("distance (mm)");
        % hold off

  
    %right insole
    
        subplot(2,2,2);
        I = imread("INSOLER.png");
        desiredRows = trialInsole.forefootWidth;
        desiredColumns = trialInsole.length;
        J = imresize(I, [desiredColumns, desiredRows]);
        %[rows, columns, numberOfColorChannels] = size(K); % double check for ts
        imshow(J)
        
        hold on
        
        j=1;
        temp = data{j}; rightInsole1 = temp{1,4};
        scatter(rightInsole1.COP_X,rightInsole1.COP_Y,colorWayNON,".");

        j = 2;
        temp = data{j}; rightInsole2 = temp{1,4};
        scatter(rightInsole2.COP_X,rightInsole2.COP_Y,colorWayNON,".");

        axis on
        title("Right Insole COP NON");
        xlabel ("distance (mm)"); ylabel ("distance (mm)");

        hold off

        % subplot(2,2,4);
        % I = imread("INSOLER.png");
        % desiredRows = trialInsole.forefootWidth;
        % desiredColumns = trialInsole.length;
        % J = imresize(I, [desiredColumns, desiredRows]);
        % %[rows, columns, numberOfColorChannels] = size(K); % double check for ts
        % imshow(J)
        % 
        % hold on
        % 
        % j = 5;
        % temp = data{j}; rightInsole = temp{1,4};
        % scatter(rightInsole.COP_X,rightInsole.COP_Y,colorWayKON,"."); alpha(a);
        % 
        % j = 6;
        % temp = data{j}; rightInsole = temp{1,4};
        % scatter(rightInsole.COP_X,rightInsole.COP_Y,colorWayKON,"."); alpha(a);
        % 
        % axis on
        % title("Right Insole COP KON");
        % xlabel ("distance (mm)"); ylabel ("distance (mm)");
        % hold off

        sgtitle("Combined SS Trials");

        % ax = axes('Position', [0.3, 0.02, 0.4, 0.08], 'Visible', 'off');
        % lgd = legend(ax, [h1, h3], {'NON', 'KON'}, 'location','southoutside');
        % lgd.NumColumns = 2;  % Display legend in two columns


% %FS Speed          
% figure(2)
%     %left insole
%         subplot(2,2,1);
%         I = imread("INSOLEL.png");
%         %imshow(I)
%         I = imrotate(I,180); %idk why the left one specifically is weird
%         desiredRows = trialInsole.forefootWidth;
%         desiredColumns = trialInsole.length;
%         K = imresize(I, [desiredColumns, desiredRows]); %scales image to right real world size
%         %[rows, columns, numberOfColorChannels] = size(K);% double check for ts
%         imshow(K)
% 
%         hold on
%         j=3;
%         temp1 = data{j}; leftInsole = temp{1,3};
%         h1 = scatter(leftInsole.COP_X, leftInsole.COP_Y,colorWayNON,".");
% 
% 
%         j = 4;
%         temp2 = data{j}; leftInsole = temp{1,3}; 
%         h2 = scatter(leftInsole.COP_X, leftInsole.COP_Y,colorWayNON,".");
% 
%         xlabel ("distance (mm)"); ylabel ("distance (mm)");
%         title("Left Insole COP NON");
%         axis on
%         hold off
% 
%         subplot(2,2,3);
%         I = imread("INSOLEL.png");
%         %imshow(I)
%         I = imrotate(I,180); %idk why the left one specifically is weird
%         desiredRows = trialInsole.forefootWidth;
%         desiredColumns = trialInsole.length;
%         K = imresize(I, [desiredColumns, desiredRows]); %scales image to right real world size
%         %[rows, columns, numberOfColorChannels] = size(K);% double check for ts
%         imshow(K)
% 
%         hold on
% 
%         j = 7;
%         temp = data{j}; leftInsole = temp{1,3}; 
% 
%         h3 = scatter(leftInsole.COP_X, leftInsole.COP_Y,colorWayKON,"."); alpha(a);
% 
%         j = 8;
%         temp = data{j}; leftInsole = temp{1,3}; 
%         h4 = scatter(leftInsole.COP_X, leftInsole.COP_Y,colorWayKON,"."); alpha(a);
% 
%         axis on;
%         title("Left Insole COP KON");
%         xlabel ("distance (mm)"); ylabel ("distance (mm)");
%         hold off
% 
% 
%     %right insole
%         subplot(2,2,2);
%         I = imread("INSOLER.png");
%         desiredRows = trialInsole.forefootWidth;
%         desiredColumns = trialInsole.length;
%         J = imresize(I, [desiredColumns, desiredRows]);
%         %[rows, columns, numberOfColorChannels] = size(K); % double check for ts
%         imshow(J)
% 
%         hold on
% 
%         j=3;
%         temp = data{j}; rightInsole = temp{1,4};
%         scatter(rightInsole.COP_X,rightInsole.COP_Y,colorWayNON,".");
% 
%         j = 4;
%         temp = data{j}; rightInsole = temp{1,4};
%         scatter(rightInsole.COP_X,rightInsole.COP_Y,colorWayNON,".");
% 
%         axis on
%         title("Right Insole COP NON");
%         xlabel ("distance (mm)"); ylabel ("distance (mm)");
%         hold off
% 
%         subplot(2,2,4);
%         I = imread("INSOLER.png");
%         desiredRows = trialInsole.forefootWidth;
%         desiredColumns = trialInsole.length;
%         J = imresize(I, [desiredColumns, desiredRows]);
%         %[rows, columns, numberOfColorChannels] = size(K); % double check for ts
%         imshow(J)
% 
%         hold on
% 
%         j = 7;
%         temp = data{j}; rightInsole = temp{1,4};
%         scatter(rightInsole.COP_X,rightInsole.COP_Y,colorWayKON,"."); alpha(a);
% 
%         j = 8;
%         temp = data{j}; rightInsole = temp{1,4};
%         scatter(rightInsole.COP_X,rightInsole.COP_Y,colorWayKON,"."); alpha(a);
% 
%         axis on
%         title("Right Insole COP KON");
%         xlabel ("distance (mm)"); ylabel ("distance (mm)");
%         hold off
% 
%     sgtitle("Combined FS Trials");
% 
%     ax = axes('Position', [0.3, 0.02, 0.4, 0.08], 'Visible', 'off');
%     lgd = legend(ax, [h1, h3], {'NON', 'KON'}, 'location','southoutside');
%     lgd.NumColumns = 2;  % Display legend in two columns
% 
% saveas(figure(1),fullfile(path_Data,'Combined'), 'jpeg');  
   
end%end function