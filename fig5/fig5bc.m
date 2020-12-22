clear;

%Here we will compute averaged time series for trajectories through the
%cases of bulk force vs invasion.  

%system size
N = 40;
M = 20;

%parameters for cell growth
h   =  5*M*N;

%time step
dt  = .0001;

%initiating propensity vector
aa  = ones(3,1);

%growth rate modulation
kappa = 0.1;

%probabilities of rotation
prob1 =  0.1;
prob2 =  0.4;

%auxiliary parameter for coarse-graining time series data
div = 1000;

%number of time points
NN = 150000;

%time vector
tt = 0:dt:NN*dt;

%number of time data points we are taking
cc = NN/div;

%number of samples
n = 1000;

%cell for storing time series data of matrices.
ZZ = cell(n,1);

%beginning of simulations
for kk = 1:n

%temp for strain fraction
Z = zeros(NN+1,1);

%states of cells at any state
numbers = [2,1,-1,-2];

%initializing matrix for lattice
A = zeros(M,N);

%Setting up A-B wall for bulk displacement
 for i = 1:N
    for j = 1:M
              if(i <= N/2)
                  A(j,i) = 2;
              else
                  A(j,i) = 1;
              end

    end
end

%Setting up random number of stripes initially for invasion
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

%this code must be kept here
count = 0;
for i = 1:M
    for k = 1:N
        if abs(A(i,k) == 2)
            count = count + 1;
        end
    end
end
        
%THE FOLLOWING CODE IS STILL USED FOR INVASION.  PLEASE UNCOMMENT!!

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

% for i = 1:N
%     for j = 1:M
%         vv = [-1, 1];
%         Q1 = randi([1,2],1);
%         Q2 = randi([1,2],1);
%         A(j,i) = vv(Q1)*Q2;
%     end
% end
% 
%CODE TILL HERE IS FOR SETUP FOR INVASION
            

%variable keeping track of fraction of smaller cells
Z(1) = count/(M*N);

%code for letting lattice evolve starts
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
   % if abs(A(R1,R2)) == 2
        if(check > 0)
            A(R1+1:M,R2) = [-A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2));A(R1+1:M-1,R2)];
        else
            A(R1,R2+1:N) = [-A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2)),A(R1,R2+1:N-1)];
        end
    %end
elseif(ind == 2)
    %if abs(A(R1,R2)) == 2
        if(check > 0)
            A(1:R1-1,R2) = [A(2:R1-1,R2); -A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2))];
        else
            A(R1,1:R2-1) = [A(R1,2:R2-1), -A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1)+(KK<prob2)*(abs(A(R1,R2))==2))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)+(KK>prob2)*(abs(A(R1,R2))==2))];
        end
    %end
end

count = 0;

for j = 1:M
    for i = 1:N
        if(abs(A(j,i)) == 2)
            count = count + 1;
        end
    end
end

Z(jjjj+1) = count/(M*N);
    
end

ZZ{kk} = Z;
end

for j = 1:n
    if ZZ{j}(end) ~= 1 
        ZZ{j}=[];
    end
end

%code to ensure all data entries are nonempty
%used to keep code clean

index = cellfun(@isempty,ZZ) == 0;
ZZZ = ZZ(index);
nN = sum(index);

AVG = zeros(NN+1,1);
for j = 1:nN
    AVG = AVG + ZZZ{j};
end
AVG = AVG/nN;



%plotting
figure(1)
plot(tt,AVG,'b-','LineWidth',4)
set(gca,'fontsize',20)
xlabel('time')
ylabel('strain fraction')
hold on

for j = 100:100:nN
    p1 = plot(tt,ZZZ{j},'b-','LineWidth',1);
    p1.Color(4) = 0.3;
end
