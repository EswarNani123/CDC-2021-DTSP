function [DubinsTourCost, alpha, ElapsedTime] = ABA(s,tour_seq_des,rho,AngleAtFirstPoint)
% %n = input('Enter number of points: '); % Inputting number of points
% iter = 100;  %input('Enter number of iterations: '); % Inputting number of iterations
alpha(1) = AngleAtFirstPoint; %input('Enter orientation angle at point 1: ');
DubinsTourCost = 0;
% rho = 10; %input('Enter minimum turning radius: ');
ss = stateSpaceDubins; % Creates an object in Dubins state space
ss.MinTurningRadius = rho; % MinTurningRadious is the propert of the object stateSpaceDubins
% n_all = 3:13; % All possible values that n takes
%for l_i = 1:length(n_all)
    n = length(s); % Taking each value of n at a time
%     txt_name = ['n' num2str(n) '_iter100_rho' num2str(rho) '.txt'];
% %     diary(txt_name)
%     s = []; % Initialization of set of points
%    for i=1:iter
%         dubins_tour_length = 0; % Initialization of dubins tour length
%         % Following code creates .tsp file
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         str1 = ['NAME: ' 'n' num2str(n)];
%         str2 = 'TYPE: ATSP'; % Change ATSP to any other thing depending on requirement. For more details, look at
%         % http://akira.ruc.dk/~keld/research/LKH/
%         str3 = ['DIMENSION: ' num2str(n)];
%         str4 = 'EDGE_WEIGHT_TYPE: EXPLICIT';
%         str5 = 'EDGE_WEIGHT_FORMAT: FULL_MATRIX' ;
%         str6 = 'EDGE_WEIGHT_SECTION';
%         xmax = 2.5; % Maximum value of x_{1}
%         xmin = -2.5; % Minimum value of x_{1}
%         ymax = 2.5; % Maximum value of x_{2}
%         ymin = -2.5; % Minimum value of x_{2}
%         s(:,1) = xmin+(xmax-xmin)*rand(n,1); % Generating the set s containing points that are distributed according to
%         s(:,2) = ymin+(ymax-ymin)*rand(n,1); % uniform probability distribution
%         [row_s,~] = size(s); % Finding number of rows and columns of s 
%         % Following code creates a cost matrix
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for z_i = 1:row_s
%             for z_j = 1:row_s
%                 if z_i ~= z_j
%                     cost(z_i,z_j) = sqrt((s(z_i,1)-s(z_j,1))^2+(s(z_i,2)-s(z_j,2))^2)*10000; % Calling the function cost_point_point
%                     % to find out the point to point time. Multiplying the
%                     % point to point time with 10000 to round the obtained time
%                     % to four decimal places
%                 end
%             end
%         end
%         cost = round(cost); % Round the number to nearest integer
%         aux_1 = max(cost,[],'all'); % Finding the maximum entry of cost matrix
%         for z_k = 1:row_s
%             cost(z_k,z_k) = aux_1+10000; % Makes the digonal elements of cost matrix too large
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % End of creating cost matrix
%         aux_2 = {str1;str2;str3;str4;str5;str6;cost}; % Create cell array using previouly defined strings and cost matrix
%         filename1 = ['n' num2str(n) '_' num2str(i) '.txt']; % Creating a file name that changes every iteration
%         filename1_tour = ['n' num2str(n) '_' num2str(i) '_tour' '.txt'];
%         fileID_1 = fopen(filename1,'w'); % Creates a text file that can be accessed using file ID fileID_1
%         fprintf(fileID_1,'%s\n',aux_2{1:6}); % Writes the data in cells 1 to 6 to a text file
%         % Following for loop writes the cost matrix to the end of text file
%         % with file ID fileID_1
%         for k = 1:n 
%             fprintf(fileID_1,'%.0f\t',aux_2{7}(k,1:n));
%             fprintf(fileID_1,'\n');
%         end
%         fclose(fileID_1); % Closes the text file created
%         % Following code creates .tsp file out of .txt file
%         filename2=strrep(filename1,'.txt','.atsp');
%         copyfile(filename1,filename2) % Copy the contents of filename1 to filename2
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % End of creating .atsp file
%         % Following code creates .par file
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         str7 = ['PROBLEM_FILE = ' filename2];
%         str8 = 'RUNS = 10';
%         str9 = ['TOUR_FILE = ' filename1_tour]; % Tour sequence is written to a text file named filename1_tour
%         fileID_2 = fopen(filename1,'w'); % Creates a text file that can be accessed using file ID fileID_2
%         fprintf(fileID_2,'%s\n',str7,str8,str9); % Writes the str7 and str8 to a text file
%         fclose(fileID_2); % Closes the text file created
%         % Following code creates .par file out of .txt file
%         filename3=strrep(filename1,'.txt','.par');
%         copyfile(filename1,filename3) % Copy the contents of filename1 to filename3
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % End of creating .par file
%         delete(filename1) % Deletes created text files which are not required
%         fileID_3 = fopen(filename1_tour,'w+');
%         system(['LKH-2' ' ' filename3]); % LKH-2 takes .par file as input and outputs 
%         % tour cost and related information on command window 
%         tour_seq_cell = textscan(fileID_3, '%d', 'MultipleDelimsAsOne',true, 'Delimiter','\n', 'HeaderLines',6);
%         % Tour sequence is extracted to tour_seq_cell from a text file with
%         % name filename1_tour
%         fclose(fileID_3);
%         tour_seq_mat = cell2mat(tour_seq_cell); % Converts from cell to matrix
%         [row_tour_seq_mat,~] = size(tour_seq_mat);
%         tour_seq_des = [tour_seq_mat(1:row_tour_seq_mat-1);1]; % Desired tour sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Following code determines the orientation angles
        tic;
        for j=1:(length(tour_seq_des)-2)
            alpha(j+1) = orientation_at_a_point(s(tour_seq_des(j),:), s(tour_seq_des(j+1),:), s(tour_seq_des(j+2),:)); 
        end
        alpha(n+1) = final_heading([s(tour_seq_des(n),:),alpha(n)],s(1,:),rho); % Orientation angle at point 1 when coming from 
        % last point in the tour sequence
        ElapsedTime = toc;
        % End of determing orientation angles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         delete(filename1_tour) % Deletes text file containg tour sequence
%         delete(filename2) % Deletes .atsp file
%         delete(filename3) % Deletes .par file
        % Followig for loop computes dubins path length from point 1 to point n
        for y_i = 1:(n-1)
            DubinsTourCost = DubinsTourCost+distance(ss,[s(tour_seq_des(y_i),:),alpha(y_i)],[s(tour_seq_des(y_i+1),:),alpha(y_i+1)]);
    %         dubConnObj = dubinsConnection;
    %         dubConnObj.MinTurningRadius = r;
    %         [pathSegObj, pathCosts] = connect(dubConnObj,[s(tour_seq_des(y_i),:),alpha(y_i)],[s(tour_seq_des(y_i+1),:),alpha(y_i+1)]);
    %         dbl = dbl+pathCosts;
    %         show(pathSegObj{1})
    %         hold on
        end
        % Dubins path length from point n to point 1 is being added to Dubins
        % path length from point 1 to point n
        DubinsTourCost = DubinsTourCost+distance(ss,[s(tour_seq_des(n),:),alpha(n)],[s(1,:),alpha(n+1)]); 
    %     dubConnObj = dubinsConnection;
    %     dubConnObj.MinTurningRadius = r;
    %     [pathSegObj, pathCosts] = connect(dubConnObj,[s(tour_seq_des(n),:),alpha(n)],[s(1,:),alpha(n+1)]);
    %     dbl = dbl+pathCosts;
    %     show(pathSegObj{1})
%         dtl(i) = dubins_tour_length; % Storing dubins tour length in each iteration
%         elap_time(i) = elap_time_1+elap_time_2; % Storing computation time in each iteration
%     end
%     if mod(l_i,2) == 1
%         subplot(2,1,1)
%         histogram(dtl,'BinWidth',1);
%         hold on
%         scatter(mean(dtl),0,'*','MarkerEdgeColor','k')
%         hold on
%     else
%         subplot(2,1,2)
%         histogram(dtl,'BinWidth',1);
%         hold on
%         scatter(mean(dtl),0,'*','MarkerEdgeColor','k')
%         hold on
%     end
%    fop = [mean(dtl) mean(elap_time) min(dtl) max(dtl)];
% %     diary off
% MeanDubinsTourCost(l_i) = mean(dtl);
% MeanElapsedTime(l_i) = mean(elap_time);
end