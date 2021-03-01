n = 13; %input('Enter the number of points: '); % User iputting number of points
iter = 100; %input('Enter number of iterations required: '); % User inputting # of iterations
for i=1:iter
    % Following code creates .atsp file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str1 = ['NAME: ' 'n' num2str(n)];
    str2 = 'TYPE: ATSP'; % Change ATSP to any other thing depending on requirement. For more details, look at
    % http://akira.ruc.dk/~keld/research/LKH/
    str3 = ['DIMENSION: ' num2str(n)];
    str4 = 'EDGE_WEIGHT_TYPE: EXPLICIT';
    str5 = 'EDGE_WEIGHT_FORMAT: FULL_MATRIX' ;
    str6 = 'EDGE_WEIGHT_SECTION';
%     xmax = 5; % Maximum value of x_{1}
%     xmin = -5; % Minimum value of x_{1}
%     ymax = 4; % Maximum value of x_{2}
%     ymin = -4; % Minimum value of x_{2}
    xmax = 2.5; % Maximum value of x_{1}
    xmin = -2.5; % Minimum value of x_{1}
    ymax = 2.5; % Maximum value of x_{2}
    ymin = -2.5; % Minimum value of x_{2}
    s(:,1) = xmin+(xmax-xmin)*rand(n,1); % Generating the set s containing points that are distributed according to
    s(:,2) = ymin+(ymax-ymin)*rand(n,1); % uniform probability distribution
    [rows,~] = size(s); % Finding number of rows and columns of s 
    % Following code creates a cost matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for z_i = 1:rows
        for z_j = 1:rows
            if z_i ~= z_j
                cost(z_i,z_j) = sqrt((s(z_i,1)-s(z_j,1))^2+(s(z_i,2)-s(z_j,2))^2)*10000; % Calling the function cost_point_point
                % to find out the point to point time. Multiplying the
                % point to point time with 10000 to round the obtained time
                % to four decimal places
            end
        end
    end
    cost = round(cost); % Round the number to nearest integer
    aux_1 = max(cost,[],'all'); % Finding the maximum entry of cost matrix
    for z_k = 1:rows
        cost(z_k,z_k) = aux_1+10000; % Makes the digonal elements of cost matrix too large
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of creating cost matrix
    aux_2 = {str1;str2;str3;str4;str5;str6;cost}; % Create cell array using previouly defined strings and cost matrix
    filename1 = ['n' num2str(n) '_' num2str(i) '.txt']; % Creating a file name that changes every iteration
    fileID_1 = fopen(filename1,'w'); % Creates a text file that can be accessed using file ID fileID_1
    fprintf(fileID_1,'%s\n',aux_2{1:6}); % Writes the data in cells 1 to 6 to a text file
    % Following for loop writes the cost matrix to the end of text file
    % with file ID fileID_1
    for k = 1:n 
        fprintf(fileID_1,'%.0f\t',aux_2{7}(k,1:n));
        fprintf(fileID_1,'\n');
    end
    fclose(fileID_1); % Closes the text file created
    % Following code creates .atsp file out of .txt file
    filename2=strrep(filename1,'.txt','.atsp');
    copyfile(filename1,filename2) % Copy the contents of filename1 to filename2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of creating .atsp file
    % Following code creates .par file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str7 = ['PROBLEM_FILE = ' filename2];
    str8 = 'RUNS = 10';
    fileID_2 = fopen(filename1,'w'); % Creates a text file that can be accessed using file ID fileID_2
    fprintf(fileID_2,'%s\n',str7,str8); % Writes the str7 and str8 to a text file
    fclose(fileID_2); % Closes the text file created
    % Following code creates .par file out of .txt file
    filename3=strrep(filename1,'.txt','.par');
    copyfile(filename1,filename3) % Copy the contents of filename1 to filename3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of creating .par file
    delete(filename1) % Deletes created text files which are not required
    tic;
    system(['LKH-2' ' ' filename3]); % LKH-2 takes .par file as input and outputs 
    % tour cost and related information on command window 
    delete(filename2)
    delete(filename3)
end