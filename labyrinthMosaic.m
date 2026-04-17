numberOfTilesInHeight = 15;

numberOfLabyrinthsToFind = 2; 
sizeOfLineInLaby = 15; 

%%

%imageFile = 'bart.png'; 
imageFile = 'lincoln.png'; 

targetImage = imread(imageFile);
[heightTarget,widthTarget,depth] = size(targetImage);

numberOfTilesInWidth = floor(widthTarget/heightTarget*(numberOfTilesInHeight+1));

heightTarget = floor(heightTarget/(numberOfTilesInHeight))*(numberOfTilesInHeight);
widthTarget  = floor(widthTarget/(numberOfTilesInWidth))*(numberOfTilesInWidth);

targetImage = imresize(targetImage,[heightTarget,widthTarget]);

pixelsInHeight = floor(heightTarget/numberOfTilesInHeight);
pixelsInWidth  = floor(widthTarget/numberOfTilesInWidth);

%% Find pixel averages over subblocks

pixelAverage = zeros(numberOfTilesInHeight,numberOfTilesInWidth);

pixelRow = 1; 
for i=1:numberOfTilesInHeight
    pixelCol = 1; 
    for j=1:numberOfTilesInWidth
        pixelAverage(i,j) = mean(mean(targetImage(pixelRow:pixelRow+pixelsInHeight-1,pixelCol:pixelCol+pixelsInWidth-1))); 
        pixelCol = pixelCol + pixelsInWidth; 
    end
    pixelRow = pixelRow + pixelsInHeight; 
end


%% Create indexing of decision variables
%  x_i,j,k is whether piece k is placed in block (i,j) where

numberOfTotalTiles = 7;
numberOfDecisionVariables = 0;
gridIndexing = zeros(numberOfTilesInHeight,numberOfTilesInWidth,numberOfTotalTiles);

indexToGrid = [];
for i=1:numberOfTilesInHeight
    for j=1:numberOfTilesInWidth
        for indexCounter=1:numberOfTotalTiles
            numberOfDecisionVariables = numberOfDecisionVariables + 1;
            gridIndexing(i,j,indexCounter) = numberOfDecisionVariables;
            indexToGrid(numberOfDecisionVariables,:) = [i,j,indexCounter];
        end
    end
end

%% Create objective function

objectiveFunction = zeros(numberOfDecisionVariables,1);
for i=1:numberOfTilesInHeight
    for j=1:numberOfTilesInWidth
        for indexCounter=1:numberOfTotalTiles-1
            index = gridIndexing(i,j,indexCounter);
            objectiveFunction(index) = pixelAverage(i,j)^2;
        end
        index = gridIndexing(i,j,numberOfTotalTiles);
        objectiveFunction(index) = (255 - pixelAverage(i,j))^2;        
    end
end

%% create constraints - one row per constraint

Aeq = [];
beq = []; constraintRowCounter = 0;

A = [];
b = []; constraintLessThanRowCounter = 0;

%% Constraint - one tile per spot

for i=1:numberOfTilesInHeight
    for j=1:numberOfTilesInWidth
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        for indexCounter=1:numberOfTotalTiles
            constraintToAdd(gridIndexing(i,j,indexCounter)) = 1;
        end
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 1;
    end
end

%% Constraint - boundaries

for i=1:numberOfTilesInHeight
    for j=1:numberOfTilesInWidth-1
        % x_i,j,1 + x_i,j,3 + x_i,j,5 = x_i,j+1,2 + x_i,j+1,4 + x_i,j+1,5
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,1)) = 1;
        constraintToAdd(gridIndexing(i,j,3)) = 1;
        constraintToAdd(gridIndexing(i,j,5)) = 1;
        constraintToAdd(gridIndexing(i,j+1,2)) = -1;
        constraintToAdd(gridIndexing(i,j+1,4)) = -1;
        constraintToAdd(gridIndexing(i,j+1,5)) = -1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
        % x_i,j,2 + x_i,j,4 + x_i,j,6 = x_i,j+1,1 + x_i,j+1,3 + x_i,j+1,6
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,2)) = 1;
        constraintToAdd(gridIndexing(i,j,4)) = 1;
        constraintToAdd(gridIndexing(i,j,6)) = 1;
        constraintToAdd(gridIndexing(i,j,7)) = 1;
        constraintToAdd(gridIndexing(i,j+1,1)) = -1;
        constraintToAdd(gridIndexing(i,j+1,3)) = -1;
        constraintToAdd(gridIndexing(i,j+1,6)) = -1;
        constraintToAdd(gridIndexing(i,j+1,7)) = -1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
    end
end

for i=1:numberOfTilesInHeight-1
    for j=1:numberOfTilesInWidth
        % x_i,j,1 + x_i,j,2 + x_i,j,6 = x_i+1,j,3 + x_i+1,j,4 + x_i+1,j,6
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,1)) = 1;
        constraintToAdd(gridIndexing(i,j,2)) = 1;
        constraintToAdd(gridIndexing(i,j,6)) = 1;
        constraintToAdd(gridIndexing(i+1,j,3)) = -1;
        constraintToAdd(gridIndexing(i+1,j,4)) = -1;
        constraintToAdd(gridIndexing(i+1,j,6)) = -1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
        % x_i,j,3 + x_i,j,4 + x_i,j,5 = x_i+1,j,1 + x_i+1,j,2 + x_i+1,j,5
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,3)) = 1;
        constraintToAdd(gridIndexing(i,j,4)) = 1;
        constraintToAdd(gridIndexing(i,j,5)) = 1;
        constraintToAdd(gridIndexing(i,j,7)) = 1;
        constraintToAdd(gridIndexing(i+1,j,1)) = -1;
        constraintToAdd(gridIndexing(i+1,j,2)) = -1;
        constraintToAdd(gridIndexing(i+1,j,5)) = -1;
        constraintToAdd(gridIndexing(i+1,j,7)) = -1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
    end
end

%% Set decision variables to zero for those that cannot appear in regions

% upper-left corner
i = 1;
j = 1;
k = 1;
constraintToAdd = zeros(1,numberOfDecisionVariables);
constraintToAdd(gridIndexing(i,j,k)) = 1;
k = 7;
constraintToAdd(gridIndexing(i,j,k)) = 1;
Aeq = [Aeq; constraintToAdd];
constraintRowCounter = constraintRowCounter + 1;
beq(constraintRowCounter) = 1;

% upper-right corner
i = 1;
j = numberOfTilesInWidth;
k = 2;
constraintToAdd = zeros(1,numberOfDecisionVariables);
constraintToAdd(gridIndexing(i,j,k)) = 1;
k = 7;
constraintToAdd(gridIndexing(i,j,k)) = 1;
Aeq = [Aeq; constraintToAdd];
constraintRowCounter = constraintRowCounter + 1;
beq(constraintRowCounter) = 1;

% lower-right corner
i = numberOfTilesInHeight;
j = numberOfTilesInWidth;
k = 4;
constraintToAdd = zeros(1,numberOfDecisionVariables);
constraintToAdd(gridIndexing(i,j,k)) = 1;
k = 7;
constraintToAdd(gridIndexing(i,j,k)) = 1;
Aeq = [Aeq; constraintToAdd];
constraintRowCounter = constraintRowCounter + 1;
beq(constraintRowCounter) = 1;

% lower-left corner
i = numberOfTilesInHeight;
j = 1;
k = 3;
constraintToAdd = zeros(1,numberOfDecisionVariables);
constraintToAdd(gridIndexing(i,j,k)) = 1;
k = 7;
constraintToAdd(gridIndexing(i,j,k)) = 1;
Aeq = [Aeq; constraintToAdd];
constraintRowCounter = constraintRowCounter + 1;
beq(constraintRowCounter) = 1;

% top-row
i = 1;
for j=1:numberOfTilesInWidth
    for k=[3 4 6]
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,k)) = 1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
    end
end

% right-row
j = numberOfTilesInWidth;
for i=1:numberOfTilesInHeight
    for k=[1 3 5]
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,k)) = 1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
    end
end

% bottom row
i = numberOfTilesInHeight;
for j=1:numberOfTilesInWidth
    for k=[1 2 6]
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,k)) = 1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
    end
end

% left-row
j = 1;
for i=1:numberOfTilesInHeight
    for k=[2 4 5]
        constraintToAdd = zeros(1,numberOfDecisionVariables);
        constraintToAdd(gridIndexing(i,j,k)) = 1;
        Aeq = [Aeq; constraintToAdd];
        constraintRowCounter = constraintRowCounter + 1;
        beq(constraintRowCounter) = 0;
    end
end

%% Set upper and lower bounds

lb = zeros(1,numberOfDecisionVariables);
ub = ones(1,numberOfDecisionVariables);
intcon = 1:length(objectiveFunction);

for numberOfLabysCreated=1:numberOfLabyrinthsToFind
    
    foundLabyrinth = 0;
    
    while ~foundLabyrinth
        
        x = intlinprog(objectiveFunction,intcon,A,b,Aeq,beq,lb,ub);
        
        % What tiles are in the current tour?
        indexOfOptArt = find(x > 10e-15);
        
        %% Establish where connections are for sub-tour detection
        %  Given (i,j) where 1 is connection to the right, -1 - left and 0 stay
        %              in the row and column, respectively. 2, indicates
        %              connections on both sides
        connect(1,:) = [ 1, 1];   % D,R
        connect(2,:) = [ 1,-1];   % D,L
        connect(3,:) = [-1, 1];   % U,R
        connect(4,:) = [-1,-1];   % U,L
        connect(5,:) = [-1, 1];   % L,R
        connect(6,:) = [-1, 1];   % U,D
        
        %% subtour detection
        allTilesFound = [];
        % The blank tiles do not create a subtour.
        allTilesFound = indexOfOptArt(find(indexToGrid(indexOfOptArt,3) == 7))';
        numberOfBlankTiles = length(allTilesFound); 
        while (length(allTilesFound) < numberOfTilesInWidth*numberOfTilesInHeight)
            remainingTiles = setdiff(indexOfOptArt,allTilesFound);
            firstTile = remainingTiles(1);
            tilesInLoop = 0;
            numberOfTilesInLoop = 0;
            nextTile = indexToGrid(firstTile(1),:);
            directionFlag = [1 -1];
            directionOfConnection = -1;
            previousStep = 1;
            while 1
                indexOfNextTile = gridIndexing(nextTile(1),nextTile(2),nextTile(3));
                nextTileType = nextTile(3);
                numberOfTilesInLoop = numberOfTilesInLoop + 1;
                tilesInLoop(numberOfTilesInLoop) = indexOfNextTile;
                
                step = zeros(1,2);
                if (nextTileType == 5)
                    nextDirection = find(directionFlag == -1*previousStep);
                    step(2) = connect(nextTileType,nextDirection);
                    previousStep = step(2);
                    directionOfConnection = 1;
                elseif (nextTileType == 6)
                    nextDirection = find(directionFlag == -1*previousStep);
                    step(1) = connect(nextTileType,nextDirection);
                    previousStep = step(1);
                    directionOfConnection = -1;
                else
                    nextDirection = find(directionFlag == directionOfConnection);
                    step(nextDirection) = connect(nextTileType,nextDirection);
                    previousStep = step(nextDirection);
                    directionOfConnection = -1*directionOfConnection;
                end
                indexOfNextTile = find(indexToGrid(indexOfOptArt,1) == nextTile(1) + step(1) & ...
                    indexToGrid(indexOfOptArt,2) == nextTile(2) + step(2));
                nextTile = indexToGrid((indexOfOptArt(indexOfNextTile,:)),:);
                if indexOfOptArt(indexOfNextTile) == firstTile
                    break;
                end
                
                %             h = plot((2*indexToGrid(tilesInLoop(numberOfTilesInLoop),2)+1)/2,-(2*indexToGrid(tilesInLoop(numberOfTilesInLoop),1)+1)/2,'*r');
                %             fprintf('%2d %2d %2d %2d\n',indexToGrid(tilesInLoop(end),:),directionOfConnection);
                %             pause
            end
            
            if (length(tilesInLoop) == numberOfTilesInWidth*numberOfTilesInHeight-numberOfBlankTiles)
                foundLabyrinth = 1;
                break;
            end
            
            %% Add a sub-loop constraint and go again.
            constraintToAdd = zeros(1,numberOfDecisionVariables);
            for k=1:length(tilesInLoop)
                constraintToAdd(tilesInLoop(k)) = 1;
            end
            A = [A; constraintToAdd];
            constraintLessThanRowCounter = constraintLessThanRowCounter + 1;
            b(constraintLessThanRowCounter) = length(tilesInLoop) - 2;
            
            allTilesFound = [allTilesFound, tilesInLoop];
        end
    end
    
    %% Draw labyrinth
    
    figure('Color',[0 0 0], 'Position', [100, 100, 1000, 1000]);
    
    for indexCounter=1:length(indexOfOptArt)
        currentGridPoint = indexToGrid(indexOfOptArt(indexCounter),:);
        i = currentGridPoint(1); j = currentGridPoint(2); k = currentGridPoint(3);
        %indexToGrid(numberOfDecisionVariables,:) = [i,j,k];
        switch k
            case 1
                xPoints = [j+1 (2*j+1)/2 (2*j+1)/2];
                yPoints = [-(2*i+1)/2 -(2*i+1)/2 -(i+1)];
            case 2
                xPoints = [j (2*j+1)/2 (2*j+1)/2];
                yPoints = [-(2*i+1)/2 -(2*i+1)/2 -(i+1)];
            case 3
                xPoints = [(2*j+1)/2 (2*j+1)/2 j+1];
                yPoints = [-i -(2*i+1)/2 -(2*i+1)/2];
            case 4
                xPoints = [(2*j+1)/2 (2*j+1)/2 j];
                yPoints = [-i -(2*i+1)/2 -(2*i+1)/2];
            case 5
                xPoints = [j j+1];
                yPoints = [-(2*i+1)/2 -(2*i+1)/2];
            case 6
                xPoints = [(2*j+1)/2 (2*j+1)/2];
                yPoints = [-i -(i+1)];
            case 7 
                continue 
        end
        plot(xPoints,yPoints,'-w','LineWidth',sizeOfLineInLaby);
        hold on
        axis off equal tight
    end
    
    currentDirectory = pwd;
    
    cd ../altmany-export_fig-4c015d5/
    eval(sprintf('export_fig ../labyrinth/results/laby%d.png',numberOfLabysCreated));
%   eval(sprintf('saveas(gcf,''../labyrinth/results/labyrinth%da.eps'')',numberOfLabysCreated));
    eval(['cd ',currentDirectory])
    
    %% Add a sub-loop constraint and go again.
    constraintToAdd = zeros(1,numberOfDecisionVariables);
    for k=1:length(tilesInLoop)
        constraintToAdd(tilesInLoop(k)) = 1;
    end
    A = [A; constraintToAdd];
    constraintLessThanRowCounter = constraintLessThanRowCounter + 1;
    b(constraintLessThanRowCounter) = length(tilesInLoop) - 2;
    
    allTilesFound = [allTilesFound, tilesInLoop];
    
    
end

