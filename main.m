%open a picture from file.
Input = imread('1.jpg');

%get the size of input image and get the size of output from use
[InputWidth, InputHeight, Depth] = size(Input);

%promptW = 'What is Width of output? ';
%OutputWidth = input(promptW);

%promptH = 'What is Height of output? ';
%OutputHeight = input(promptH);

%Should not to gray


%crop image(256=>400?)->
%store subimage as a 4-D array->

%preallocation

t = int32(InputWidth/4);
re = mod(t, 6);
if(re ~= 0)
    t = t+(6-re);
end

jumpPixel = 5;
Ssize = ((InputWidth-t)/jumpPixel)*((InputHeight-t)/jumpPixel); 
S = zeros(t, t, 3, Ssize);

original = zeros(512, 512, 3);
for i=1:2
    for j=1:2
        original(1+256*(i-1):256*i, 1+256*(j-1):256*j, :) = Input;
    end
end

for i = 1:((InputWidth-t)/jumpPixel)
    for j = 1:((InputHeight-t)/jumpPixel)
        S(:, :, :, (i-1)*((InputWidth-t)/jumpPixel)+j) = imcrop(Input, [1+jumpPixel*i 1+jumpPixel*j t-1 t-1]);
    end
end

%neighbor it by SSD->
%still have best path to determine

numOfoutputPart = 10;
w = ((t*5)/6) * numOfoutputPart + t/6;
h = ((t*5)/6) * numOfoutputPart + t/6;
length = t*5/6;
overlapLength = t/6;
preGoal = zeros(w, h, 3);
Goal = zeros(w, h, 3);

Store = zeros(numOfoutputPart.^2);
for i = 1:numOfoutputPart %trace the whole subimage set
    for j = 1:numOfoutputPart
        if(i==1) && (j==1) %initial goal image
            preGoal(1:t, 1:t, :) = S(:, :, :, 1);
            Goal(1:t, 1:t, :) = S(:, :, :, 1);
            Store(1) = 1;
        else
            SSD = [0 0]; %[index, SSD]
            LeftSSD = 0;
            UpSSD = 0;
            for k = 1:(InputWidth-t)/jumpPixel %Search for best SSD
                for l = 1:(InputHeight-t)/jumpPixel
                    if(i~=1) %need to check UpSSD
                        temp1 = S(length:length+overlapLength-1, 1:t, :, Store((i-1-1)*4+j));
                        temp2 = S(1:1+overlapLength-1, 1:t, :, ((k-1)*(InputHeight-t)/jumpPixel)+l) ;
                        X = temp2 - temp1;
                        UpSSD = sum( X(:).^2);
                    end
                    
                    if(j~=1) %need to check LeftSSD
                        temp1 = S(1:t, length:length+overlapLength-1, :, Store((i-1)*4+j-1));
                        temp2 = S(1:t, 1:1+overlapLength-1, :, ((k-1)*(InputHeight-t)/jumpPixel+l) );
                        Y = temp2 - temp1;
                        LeftSSD = sum( Y(:).^2);
                    end
                    
                    if( (UpSSD+LeftSSD)<SSD(2) || SSD(1) == 0) %Find best SSD
                        SSD(1) = (k-1)*((InputWidth-t)/jumpPixel)+l;
                        SSD(2) = UpSSD+LeftSSD;
                        Store( (i-1)*numOfoutputPart+j ) = SSD(1);
                    end
                end
            end
            %finish finding best matching blocks
            
            %only neighbor on SSD
            
            preGoal(1+(i-1)*length:t+(i-1)*length, 1+(j-1)*length:t+(j-1)*length, :) = S(:, :, :, Store((i-1)*numOfoutputPart+j) );
            
            %Neighboring
            errorMatrixH = zeros(overlapLength, t);
            adjacentErrorMatrixH = zeros(overlapLength, t);
            bestPathH = zeros(2, t);
            if(i ~= 1) %need to check upper overlap
                %determine the overlap pointer first (vertical)
                S1 = Goal(length*(i-2)+1:length*(i-2)+t, length*(j-1)+1:length*(j-1)+t, :);
                S2 = S(:, :, :, Store((i-1)*4+j));
                
                
                %find e(i,j) & E(i,j)
                tempH = S1(length+1:t, :, :) - S2(1:overlapLength, :, :);
                for k = 1:overlapLength
                    for l = 1:t
                        errorMatrixH(k, l) = sum(tempH(k, l).^2);
                    end
                end
                
                %initial the E matrix
                adjacentErrorMatrixH(1, :) = errorMatrixH(1, :);
                
                %store minimum error
                
                bestPathH(1, :) = adjacentErrorMatrixH(1, :);
                for k = 1:t
                bestPathH(2, k) = 1;
                end
                for k = 2:overlapLength
                    for l = 1:t
                        Min = 0;
                        if(l == 1)
                            Min = min(adjacentErrorMatrixH(k-1, l), adjacentErrorMatrixH(k-1, l+1));
                        elseif(l == t)
                            Min = min(adjacentErrorMatrixH(k-1, l-1), adjacentErrorMatrixH(k-1, l));
                        else
                            Min = min([adjacentErrorMatrixH(k-1, l-1), adjacentErrorMatrixH(k-1, l), adjacentErrorMatrixH(k, j+1)]);
                        end
                        
                        adjacentErrorMatrixH(k, l) = errorMatrixH(k, l) + Min;
                        
                        if(adjacentErrorMatrixH(k, l) < bestPathH(1, l))
                            bestPathH(1, l) = adjacentErrorMatrixH(k, l);
                            bestPathH(2, l) = k;
                        end
                    end
                end
                
            end
            
            errorMatrixW = zeros(t, overlapLength);
            adjacentErrorMatrixW = zeros(t, overlapLength);
            bestPathW = zeros(t, 2);
            if(j ~= 1)%need to check left overlap
                
                %determine the overlap pointer first (vertical)
                S1 = Goal(length*(i-1)+1:length*(i-1)+t, length*(j-2)+1:length*(j-2)+t, :);
                S2 = S(:, :, :, Store((i-1)*4+j));
                
                
                %find e(i,j) & E(i,j)
                tempW = (S1(:, 1+length:t, :) - S2(:, 1:overlapLength, :));
                for k = 1:t
                    for l = 1:overlapLength
                        errorMatrixW(k, l) = sum(tempW(k, l).^2);
                    end
                end
                
                
                %initial the E matrix
                adjacentErrorMatrixW(:, 1) = errorMatrixW(:, 1);
                
                %store minimum error
                
                bestPathW(:, 1) = adjacentErrorMatrixW(:, 1);
                for k = 1:t
                bestPathW(k, 2) = 1;
                end
                for l = 2:overlapLength
                    for k = 1:t
                        Min = 0;
                        if(k == 1)
                            Min = min(adjacentErrorMatrixW(k, l-1), adjacentErrorMatrixW(k+1, l-1));
                        elseif(k == t)
                            Min = min(adjacentErrorMatrixW(k-1, l-1), adjacentErrorMatrixW(k, l-1));
                        else
                            Min = min([adjacentErrorMatrixW(k-1, l-1), adjacentErrorMatrixW(k, l-1), adjacentErrorMatrixW(k+1, l-1)]);
                            
                        end
                        
                        adjacentErrorMatrixW(k, l) = errorMatrixW(k, l)+ Min; 
                        
                        if(adjacentErrorMatrixW(k, l) < bestPathW(k, 1))
                            bestPathW(k, 1) = adjacentErrorMatrixW(k, l);
                            bestPathW(k, 2) = l;
                        end
                    end
                end
                
            end
            
            for k = 1:t
                for l = 1:t
                    if(k > bestPathH(2, l) && l > bestPathW(k, 2) )
                        Goal( (i-1)*length+k, (j-1)*length+l, :) = S2(k , l, :);
                    end
                end
            end
        end
    end
end

%preGoal = preGoal/255;
%image(preGoal);
%imshow(preGoal);


Goal = Goal/255;
image(Goal);
imshow(Goal);


%GoalFinal = GoalFinal/255;
%image(GoalFinal);
%imshow(GoalFinal);

















