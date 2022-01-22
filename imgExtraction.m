function [tempK, tempF, time] = imgExtraction(fileName)
% Function to automatically extract data from graph images of temperature
% variation through different tile locations graphs given by
% NASA's Aeroheating Flight Experiment
% (https://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-aiaa-2001-0352.pdf)

% W Powell  06/04/21
% Required input arguments:
% fileName   - image name of tile graph at desired location

% Output arguments:
% tempK    - temperature data in Kelvin
% tempF    - temperature data in Farenheit
% time     - time data

% For example, to perform function for tile 597, with file name 597.png:
% [tempK, tempF, time] = imgExtraction('597');

% threshold to detect red/black in image
redThreshold = 80;
blackThreshold = 20; 

% initial conditions that are no available from image plot
initTempK = 292;
initTempF = 66;
initTime = 0;

img = imread(['ShuttleImgs/' fileName '.png']); % load image from file

imgNoPlot = 255 - round((img(:,:,2) + img(:,:,3))); % removes red from image
imgAxes = imgNoPlot > blackThreshold; % convert to binary image of Axes

% xAxis and yAxis locator:

yAxis = 0; % yAxis is the column of img where the y axis is located
maxSumColumn = 0; % arbitary high value to ensure sumColumn is greater on first column
sumRow = zeros(size(img,1), 1);
prevCell = 1;
yLimCoord = [];

% loops through each cell and detects x axis and y axis by summing pixel
% values for column and row respectively
for column = 1:size(imgAxes,2)
    sumColumn = 0; % reset sum of pixels for new column to zero
    for row = 1:size(imgAxes,1)
        cell = double(imgAxes(row,column)); % finds value of current cell
        if prevCell == 0 && cell == 1 % store all cells that transition from white to black
            yLimCoord(end+1,:) = [row column];
        end
        sumRow(row,1) = sumRow(row,1) + cell;
        sumColumn = sumColumn + cell; % sums all cell values
        prevCell = cell;
    end
    if sumColumn > maxSumColumn % axes are black therefore lowest sum of pixels will be at y axis column
        yAxis = column;
        maxSumColumn = sumColumn; % updates new largest summed column
    end
    
end

[~, xAxis] = max(sumRow); % sets xAxis to maximum pixel sum for the row (when row is black due to axes)
prevCell = 0;
xLim = 0;
% locates x axis limit
for column = 1:size(imgAxes,2)
    cell = imgAxes(xAxis,column);
    if prevCell == 1 && cell == 0 % store all cells that transition from white to black
            xLim = column - 1;
    end
    prevCell = cell;
end
yLim = yLimCoord(find((yLimCoord(:,2) == yAxis),1),1); % locates y limit

% used for graph transformations
setAxesTemp = [yLim, xAxis]; 
setAxesTime = [xLim, yAxis];

% crop image to axes
imgCrop = img((yLim : xAxis),(yAxis : xLim),:);
% extracts red plot line from image
imgRed = imgCrop(:,:,1) - round((imgCrop(:,:,2) + imgCrop(:,:,3)));
% converts image to logic matrix / binary image 
imgBinary = imgRed > redThreshold; 
imgBinary = flip(imgBinary); % flips matrix

avgYCoord = zeros(1,size(imgBinary,2)); % pre allocates space for vector
% averages all y values in each column, since multiple pixels in each column
for i = 1:size(imgBinary,2)
    sumYCoord = 0;
    sumColumn = 0;
    for j = 1:size(imgBinary,1)
        if imgBinary(j,i) == 1
                sumColumn = sumColumn + 1;
                sumYCoord = sumYCoord + j;
                avgYCoord(1,i) = (sumYCoord / sumColumn);
        end
    end
end

sizeYCoord = size(avgYCoord);
x = 1:sizeYCoord(2); % x axis of image is number of columns

% graph transformations
tempF = ((2000/(xAxis - yLim)) .* (avgYCoord)); % temp data in fahrenheit
time = ((2000 / (xLim-yAxis))  .* (x)); % time data in seconds

tempK = (tempF - 32) .* (5/9)+ 273.15; % temp data in kelvin

% initial conditions which are unkown from image
tempK(1) = initTempK;
tempF(1) = initTempF;
tempK(avgYCoord == 0) = initTempK; % set unkown data points at low time values
tempF(avgYCoord == 0) = initTempF; 
time(1) = initTime;