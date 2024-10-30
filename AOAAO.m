
function [Best_FF,Best_P,Conv_curve]=AOAAO(N,M_Iter,LB,UB,Dim,F_obj)
display('AOAAO Working');
%Two variables to keep the positions and the fitness value of the best-obtained solution

Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);

%Initialize the positions of solution
X=initialization(N,Dim,UB,LB);
Xnew=X;
Ffun=zeros(1,size(X,1));% (fitness values)
Ffun_new=zeros(1,size(Xnew,1));% (fitness values)

MOP_Max=1;
MOP_Min=0.2;
C_Iter=1;
Alpha=5;
Mu=0.499;
chaos = PiecewiseLinear(M_Iter);

for i=1:size(X,1)
    Ffun(1,i)=F_obj(X(i,:));  
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
    to = 1:Dim;
    u = .0265;
    r0 = 10;
    r = r0 +u*to;
    omega = .005;
    phi0 = 3*pi/2;
    phi = -omega*to+phi0;
    x = r .* sin(phi);  % Eq. (9)
    y = r .* cos(phi); % Eq. (10)
    QF=C_Iter^((2*rand()-1)/(1-M_Iter)^2); % Eq. (15)  
    

while C_Iter<M_Iter+1  %Main loop
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio 
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function 
    E1=2*(1-(C_Iter/M_Iter)); 
    %Update the Position of solutions
    for i=1:size(X,1)   % if each of the UB and LB has a just value 
        E0=2*rand()-1; 
        Escaping_Energy=E1*(E0)+(0.9-0.4*C_Iter/M_Iter)*chaos(C_Iter);  
        
        if abs(Escaping_Energy)>=1
            
 
            q=rand(); 
            if q<0.5
                X1=Best_P(1,:)*(1-C_Iter/M_Iter)+(mean(X(i,:))-Best_P(1,:))*rand(); % Eq. (3) and Eq. (4)
                if F_obj(X1)<F_obj(X(i,:))
                    Xnew(i,:)=X1;
                end
            elseif q>=0.5
                X2=Best_P(1,:).*Levy(Dim)+X((floor(N*rand()+1)),:)+(y-x)*rand;       % Eq. (5)
                if F_obj(X2)<F_obj(X(i,:))
                    Xnew(i,:)=X2;
                end
            end
            
        elseif abs(Escaping_Energy)<1
           r1=rand();
            if (size(LB,2)==1)
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        X3=Best_P(1,:)/(MOP+eps)*((UB-LB)*Mu+LB);
                       if F_obj(X3)<F_obj(X(i,:))
                        Xnew(i,:)=X3;
                        end
                    else
                        X4=Best_P(1,:)*MOP*((UB-LB)*Mu+LB);
                       if F_obj(X4)<F_obj(X(i,:))
                        Xnew(i,:)=X4;
                        end
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,:)=Best_P(1,:)-MOP*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,:)=Best_P(1,:)+MOP*((UB-LB)*Mu+LB);
                    end
                end               
            end    
            
        end
        
        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
 
        Ffun_new(1,i)=F_obj(Xnew(i,:)); % calculate Fitness function 
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
        end
       
    end
    

    %Update the convergence curve
    Conv_curve(C_Iter)=Best_FF;
    
    %Print the best solution details after every 50 iterations
    if mod(C_Iter,50)==0
        display(['At iteration ', num2str(C_Iter), ' the best solution fitness is ', num2str(Best_FF)]);
    end
     
    C_Iter=C_Iter+1;  % incremental iteration
   
end
end
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end

function [xlist] =PiecewiseLinear(maxIter)
    initial_value=0.7;
    a=0.6;
    x = initial_value;
    xlist = zeros(1,maxIter);
    for i=1:maxIter
        if x > 0 && x < 1-a
                x = x/(1-a);
           
        else
                x = (x-(1-a))/a;
         end
        xlist(i) = x;
    end
end



