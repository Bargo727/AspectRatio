clear;

%In this program we illustrate the existence of invasion and bulk
%displacement in the lattice model (LM)

%lattice size

%columns
N = 40;
%rows
M = 20;

%parameters for simulation

%growth rate of cells
h   =  5*M*N;

%time step
dt  = .0001;

%propensity vector
aa  = ones(3,1);

%parameter capturing mechanical forces effect on cell growth rate
kappa = 0.1;

%probability of rotation for a newly born cell of strain 1 or 2
prob1 = .01;
prob2 = .5;

%parameter for establishing at what time snaps the full lattice is output
%as a figure
div = 1000;

%number of time steps
NN = 50000;

%setting up cells to store lattice microscopic configurations to store and
%produce as a movie if desired
cc = NN/div;
EE = cell(cc+1,1);

%vector storing states of cells in lattice:
%2:  cell of strain 2 and vertically oriented
%-2: cell of strain 2 and horizontally oriented
%1:  cell of strain 1 and vertically oriented
%-1: cell of strain 1 and horizontally oriented

numbers = [2,1,-1,-2];

%matrix that represents lattice
A = zeros(M,N);

%setting up A-B wall for bulk displacement illustration
for i = 1:N
    for j = 1:M
              if(i <= N/2)
                  A(j,i) = 2;
              else
                  A(j,i) = 1;
              end

    end
end

%setting up random strip distribution for invasion illustration
% R = randi([1,N-1],randi([1,N],1,1),1);
% R = unique(R);
% LR = length(R);
% mR = min(R);
% MR = max(R);
% 
% for j = 1:LR
%     if j == 1
%         A(:,1:R(j)) = randi([1,2],1,1)*ones(M,R(j));
%     else
%         if A(M/2,R(j-1)) == 1
%             A(:,R(j-1)+1:R(j)) = 2*ones(M,R(j) - R(j-1));
%         else
%             A(:,R(j-1)+1:R(j)) = ones(M,R(j) - R(j-1));
%         end
%     end
% end
% 
% if A(M/2,R(end)) == 1
%     A(:,R(end)+1:N) = 2*ones(M,N-R(end));
% else
%     A(:,R(end)+1:N) = ones(M,N-R(end));
% end
% 
%            
% count = 0;
% for i = 1:M
%     for k = 1:N
%         if abs(A(i,k) == 2)
%             count = count + 1;
%         end
%     end
% end
%         
% 
% while count ~= M*N/2
%     if count > M*N/2
%         ii = find(sum(A) == 2*M, 1);
%         A(:,ii) = ones(M,1);
%         
%         count = 0;
%         for i = 1:M
%             for k = 1:N
%                 if abs(A(i,k) == 2)
%                     count = count + 1;
%                 end
%             end
%         end
% 
%     elseif count < M*N/2
%         ii = find(sum(A) == M,1);
%         A(:,ii) = 2*ones(M,1);
%         
%         count = 0;
%         for i = 1:M
%             for k = 1:N
%                 if abs(A(i,k) == 2)
%                     count = count + 1;
%                 end
%             end
%         end
%     end
% end
%
%all code up till here codes for random stripe distribution

%first cell displays initial conditions
EE{1} = A;

%code to output the initial data
figure(1)
        for i = 1:N
        for j=1:M
            temp = A(j,i);
            thisStrain = abs(temp);
            if thisStrain==1 
                cellPaint = 'b'; 
            else
                cellPaint = 'g'; 
            end
            x = 3*i;%jsonData.frames(i).cells(j).d(1); xi = round(x) + 1;
            y = 3*j;%jsonData.frames(i).cells(j).d(2); yj = round(y) + 1;
            phi = (pi/2)*(temp >0) + pi*(temp<0);
            l = 3;%jsonData.frames(i).cells(j).d(4);
            %comp = jsonData.frames(i).cells(j).d(5);
            cpoly = makeCell(x,y, phi,l);
            fill(cpoly(1,:), cpoly(2,:), cellPaint); hold on;    
        end
        end
            hold off;
            set(gca,'fontsize',20)
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            axis([-1 3*N+1 -1 3*M+1]);
            axis equal;
         


%loop to begin simulation
for jjjj = 1:NN
    R1 = randi([1,M]);
    R2 = randi([1,N]);
    check = A(R1,R2);
    aa(1) = (check>0)*h*(R1~=M)*exp(-kappa*(M-R1))*dt + (check<0)*h*(R2~=N)*exp(-kappa*(N-R2))*dt;
    aa(2) = (check>0)*h*(R1~=1)*exp(-kappa*(R1-1))*dt + (check<0)*h*(R2~=1)*exp(-kappa*(R2-1))*dt;
    aa(3) = 1 - aa(1) - aa(2);
    
    selection_intv1 = cumsum(aa); 
    selection_intv1 = selection_intv1./selection_intv1(end);
    ind = find(selection_intv1 > rand,1);
    
    

KK = rand;
if(ind == 1)
        if(check > 0)
            A(R1+1:M,R2) = [-A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2));A(R1+1:M-1,R2)];
        else
            A(R1,R2+1:N) = [-A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2)),A(R1,R2+1:N-1)];
        end
elseif(ind == 2)
        if(check > 0)
            A(1:R1-1,R2) = [A(2:R1-1,R2); -A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2))];
        else
            A(R1,1:R2-1) = [A(R1,2:R2-1), -A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2))];
        end
end

    if mod(jjjj,div)==0
        EE{(jjjj/div)+1} = A;
    end
    
    if mod(jjjj,1000) == 0
        figure(jjjj)
        for i = 1:N
        for j=1:M
            temp = A(j,i);
            thisStrain = abs(temp);
            if thisStrain==1 
                cellPaint = 'b'; 
            else
                cellPaint = 'g'; 
            end
            x = 3*i;%jsonData.frames(i).cells(j).d(1); xi = round(x) + 1;
            y = 3*j;%jsonData.frames(i).cells(j).d(2); yj = round(y) + 1;
            phi = (pi/2)*(temp >0) + pi*(temp<0);
            l = 3;%jsonData.frames(i).cells(j).d(4);
            %comp = jsonData.frames(i).cells(j).d(5);
            cpoly = makeCell(x,y, phi,l);
            fill(cpoly(1,:), cpoly(2,:), cellPaint); hold on;    
        end
        end
            hold off;
            set(gca,'fontsize',20)
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            axis([-1 3*N+1 -1 3*M+1]);
            axis equal;
    end
    
end

%uncomment to make video

%mf = figure(2);  clear vidfile; 
% for k=1:cc+1
%     clf;
%     for i = 1:N
%         for j=1:M
%             temp = EE{k}(j,i);
%             thisStrain = abs(temp);
%             if thisStrain==1 
%                 cellPaint = 'b'; 
%             else
%                 cellPaint = 'g'; 
%             end
%             x = 3*i;%jsonData.frames(i).cells(j).d(1); xi = round(x) + 1;
%             y = 3*j;%jsonData.frames(i).cells(j).d(2); yj = round(y) + 1;
%             phi = (pi/2)*(temp >0) + pi*(temp<0);
%             l = 3;%jsonData.frames(i).cells(j).d(4);
%             %comp = jsonData.frames(i).cells(j).d(5);
%             cpoly = makeCell(x,y, phi,l);
%             fill(cpoly(1,:), cpoly(2,:), cellPaint); hold on;    
%         end
%      end
%     hold off;
%     set(gca,'fontsize',20)
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%      axis([-1 3*N+1 -1 3*M+1]);
%      axis equal;
%      drawnow;  %will slow down the code
%     F(k) = getframe;
% end
% 
% vidfile = VideoWriter('illustration1.mp4','MPEG-4');
% vidfile.FrameRate = 2.5;
% open(vidfile);
% writeVideo(vidfile,F)
% close(vidfile)
% % 
% 
% 

