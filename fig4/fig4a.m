clear;

%In this program we investigate how p_rot changes the bulk disorder
%paramater, characterized by
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               q = (<cos(2\phi)>^2 + <sin(2\phi)>^2)^(1/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%where \phi is the angle the cells make with the horizontal.  We only
%consider isogenic populations, akin to the ABM.

%lattice size

%columns
N = 80;
%rows
M = 40;

%parameters for simulation

%growth rate of cells
h   =  5*M*N;

%time step
dt  = .0001;

%number of time steps
NN = 200000;

%parameter for establishing at what time snaps the full lattice is output
%as a figure
div = 1;

%propensity vector
aa  = ones(3,1);

%parameter capturing mechanical forces effect on cell growth rate
kappa = 0.5;

P1 = 0.0:0.05:1;
SZ = length(P1);
Q = zeros(SZ,1);
V = zeros(SZ,1);

%probability of rotation for a newly born cell of strain 
for bb = 1:SZ

prob1 = P1(bb);

%size of storage vector
sz = NN/2;

%vector storing q over time
%q_vector = zeros(sz,1);
q_vector = [];


%vector storing states of cells in lattice:
%1:  vertically oriented cell
%-1: horizontally oriented cell

numbers = [1,-1];

%matrix that represents lattice
A = zeros(M,N);


%the following code initializes the lattice to have a columnar initial
%condition
for i = 1:N
    for j = 1:M
        A(j,i) = 1;
    end
end

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
            A(R1+1:M,R2) = [-A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1));A(R1+1:M-1,R2)];
        else
            A(R1,R2+1:N) = [-A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1)),A(R1,R2+1:N-1)];
        end
elseif(ind == 2)
        if(check > 0)
            A(1:R1-1,R2) = [A(2:R1-1,R2); -A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1))];
        else
            A(R1,1:R2-1) = [A(R1,2:R2-1), -A(R1,R2)*((KK<prob1)*(abs(A(R1,R2))==1))+A(R1,R2)*((KK>prob1)*(abs(A(R1,R2))==1))];
        end
end

if (jjjj > NN/2) == 0
    sumc = 0;
    sums = 0;
    for i = 1:N
        for j = 1:M
            if A(j,i) > 0
                sumc = sumc + cos(pi);
                sums = sums + sin(pi);
            elseif A(j,i) < 0
                sumc = sumc + cos(0);
                sums = sums + sin(0);
            end
        end
    end
    sumc = sumc/(M*N);
    sums = sums/(M*N);
    q_vector = [q_vector;sqrt(sumc^2 + sums^2)]; 
end
    
    
end
Q(bb) = mean(abs(q_vector));
V(bb) = std(abs(q_vector));
end

figure(2)
errorbar(Q,P1,V,'b-o','LineWidth',3,'MarkerFaceColor','b','MarkerSize',10)
xlabel('p_{rot}')
ylabel('q')
set(gca,'fontsize',20)
