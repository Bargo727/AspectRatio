clear;

%In this code we show spatial oscillation of the interface in an AB wall

%spatial dimensions

N = 30;
M = 15;

%parameters for simulation

h      = 10*M*N;
kappa  = 0;
eps1   = 0.1;
eps2   = 0.1;
H      = 2;

%total number of rxns

NN = 10000;
div = 100;
cc = NN/div;

%time step
dt  = .0001;



%setting up matrices
%A contains the state of the cell occupying the ijth site .


A = zeros(M,N);
for i = 1:N
    for j = 1:M
        if(i<= N/2)
            A(j,i) = 2;
        else
            A(j,i) = 1;
        end
    end
end

count1 = M*N/2;
count2 = M*N/2;


prob1 = 0.01;
prob2 = 0.2;
prov1 = eps1;
prov2 = eps2;

figure(1)
        for i = 1:N
        for j=1:M
            temp = A(j,i);
            thisStrain = abs(temp);
            if thisStrain==1 
                cellPaint = [0, 114, 178]/256; 
            else
                cellPaint = [213, 94, 0]/256; 
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


%Vector to store propensities.  Because there are two possible growth
%directions for each cell, we must store them as separate propensities in a
%given vector.  Within this we must also account for the different types of
%cells, namely of strain 1 or strain 2.  Therefore, the size
%will be M x 2N.  

p = zeros(M,2*N);

for k = 1:NN
    if(count1 <= M*N/3)
        prob1 = eps1 + (1 - eps1)*((count2.^H)./(1 + count2.^H));
        prob2 = 0;%+ (1 - eps2)*((count1.^H)./(1 + count1.^H));
        prov2 = 1;
        prov1 = 1-prob1;
        
%         prob1 = 0;
%         prob2 = 0;
%         prov1 = 1;
%         prov2 = 1;

    elseif(count2 <= M*N/3)
        prob1 = 0;%eps1 + (1 - eps1)*((count2.^H)./(1 + count2.^H));
        prob2 = eps2 + (1 - eps2)*((count1.^H)./(1 + count1.^H));
        prov2 = 1-prob2;
        prov1 = 1;

    end

for j = 1:M
    for i = 1:2
        if(i == 1)
            p(j,i) = ((A(j,1) == 1)+ (A(j,1) == 2))*(j ~= M)*h*exp(-kappa*(M-j)) + ((A(j,1) == -1)+(A(j,1) == -2))*h*exp(-kappa*(N-1)); 
        elseif(i == 2)
            p(j,i) = ((A(j,1) == 1)+(A(j,1) == 2))*(j ~= 1)*h*exp(-kappa*(j-1));
        end
    end
end

for j = 1:M
    for i = 3:2*N
        if(mod(i,2)==1)
            p(j,i) = ((A(j,(i-1)/2) == 1)+ (A(j,(i-1)/2) == 2))*(j ~= M)*h*exp(-kappa*(M-j)) + ((A(j,(i-1)/2) == -1)+(A(j,(i-1)/2) == -2))*(((i-1)/2)~=N)*h*exp(-kappa*(N-((i-1)/2))); 
        elseif(mod(i,2) == 0)
            p(j,i) = ((A(j,(i-2)/2) == 1)+(A(j,(i-2)/2) == 2))*(j ~= 1)*h*exp(-kappa*(j-1)) + ((A(j,(i-2)/2) == -1)+(A(j,(i-2)/2) == -2))*(((i-2)/2)~=1)*h*exp(-kappa*(((i-2)/2)-1));
        end
    end
end


P = reshape(p',2*M*N,1);

%generating a time step for a reaction to occur based on exponentially
%distributed time

% delta_t = -log(rand(1,1))/sum(P);
% 
% T(k+1) = T(k) + delta_t;

%Selecting valid rxns only
valid_inds = P>0;
valid_P    = P(valid_inds);

%Normalizing for reaction selection process

selection_interval = cumsum(valid_P);
selection_interval = selection_interval/selection_interval(end);

%Selecting index of reaction that will occur
PP = cumsum(P)./sum(P);
tol = 0.0005;
selected_ind1 = find(selection_interval > rand(1,1),1);
RXN = selection_interval(selected_ind1);
selected_ind  = find(abs(PP-RXN) < tol,1);

%Finding index of location corresponding to the reaction
location_index = ceil(selected_ind/2);

%Finding location of cell that will induce change
r = ceil(location_index/N);
c = mod(location_index,N)*(mod(location_index,N)~=0) + N*(mod(location_index,N) ==0);

%Indexing reaction to know what type of reaction it is
rxn_ind = mod(selected_ind,2);
% if (rxn_ind == 8)
%     stop
% end

if(rxn_ind == 1)
    if A(r,c) > 0;
        R = rand(1,1);
        A(r+1:end,c) = [A(r,c)*((abs(A(r,c))==2)*(prob2 < R) + (abs(A(r,c))==1)*(prob1 < R)) - A(r,c)*((abs(A(r,c))==2)*(prob2 > R) + (abs(A(r,c))==1)*(prob1 > R));A(r+1:end-1,c)];
    else
        R = rand(1,1);
        A(r,c+1:end) = [A(r,c)*((abs(A(r,c))==2)*(prov2 < R) + (abs(A(r,c))==1)*(prov1 < R)) - A(r,c)*((abs(A(r,c))==2)*(prov2 > R) + (abs(A(r,c))==1)*(prov1 > R)),A(r,c+1:end-1)];
    end
elseif(rxn_ind == 0)
     if A(r,c) > 0;
        R = rand(1,1);
        A(1:r-1,c) = [A(2:r-1,c);A(r,c)*((abs(A(r,c))==2)*(prob2 < R) + (abs(A(r,c))==1)*(prob1 < R)) - A(r,c)*((abs(A(r,c))==2)*(prob2 > R) + (abs(A(r,c))==1)*(prob1 > R))];
     else
        R = rand(1,1);
        A(r,1:c-1) = [A(r,2:c-1),A(r,c)*((abs(A(r,c))==2)*(prov2 < R) + (abs(A(r,c))==1)*(prov1 < R)) - A(r,c)*((abs(A(r,c))==2)*(prov2 > R) + (abs(A(r,c))==1)*(prov1 > R))];
     end
end

for j = 2:M-1
    for i = 2:N-1
        if abs(A(j,i)) == 1
            c1 = 0;
            if( abs(A(j+1,i)) ~= 1)
                c1 = c1+1;
            end
            if( abs(A(j-1,i)) ~= 1)
                c1 = c1+1;
            end
            if( abs(A(j,i+1)) ~= 1)
                c1 = c1+1;
            end
            if( abs(A(j,i-1)) ~= 1)
                c1 = c1+1;
            end
            if(c1 >= 3)
                A(j,i) = sign(A(j,i))*2;
            end
        else
            c1 = 0;
            if( abs(A(j+1,i)) ~= 2)
                c1 = c1+1;
            end
            if( abs(A(j-1,i)) ~= 2)
                c1 = c1+1;
            end
            if( abs(A(j,i+1)) ~= 2)
                c1 = c1+1;
            end
            if( abs(A(j,i-1)) ~= 2)
                c1 = c1+1;
            end
            if(c1 >= 3)
                A(j,i) = sign(A(j,i));
            end
        end
    end
end


count1 = 0;
count2 = 0;

for j = 1:M
    for i = 1:N
        if (abs(A(j,i))== 1)
            count1 = count1 + 1;
        elseif (abs(A(j,i)) == 2)
            count2 = count2 + 1;
        end
    end
end

 if mod(k,500) == 0
        figure(k)
        for i = 1:N
        for j=1:M
            temp = A(j,i);
            thisStrain = abs(temp);
            if thisStrain==1 
                cellPaint = [0, 114, 178]/256; 
            else
                cellPaint = [213, 94, 0]/256; 
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

