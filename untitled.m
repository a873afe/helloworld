clear all
%%% parameter setting 
    d0     = 0.7      ;   
    Mc     = 10:20:200;     % # of UEs       
    N      = 64       ;     % # of Preambles
    K      = 16       ;     % # of BS antennas
    
    n_pw   = 0.01     ;     % noise power is arbitrarily set
    Power1 = 0.3      ;
    Power2 = 1 - Power1;

    Msize  = length(Mc);
    
    noma_c = zeros(1,Msize);
     
    oma_c  = zeros(1,Msize);
    
for m = 1:Msize
       
    M= Mc(m);
    simulation_times = 20000;

    P1res = zeros(1,simulation_times);
    P2res = zeros(1,simulation_times);
    P3res = zeros(1,simulation_times);  
  
    for j = 1 : simulation_times 
       
       % fprintf('loop = %f %%\n', 100*j/simulation_times)
        %% UE map building 
        UEmap = zeros(M,3); 
        for i = 1:M
            r = sqrt(rand);
            theta = rand*(2*pi);
            x = r * cos(theta);
            y = r * sin(theta);
            UEmap(i,:) = [x  y  sqrt(power(x,2)+power(y,2))];
        end
    %UEmap
    
    %% Random Channel 
   
    H_temp = Rayleigh_Fading(K,M);
   
    H = zeros(K,M);
    for i = 1:M % loss
        H(:,i) = H_temp(:,i)*(sqrt(1-UEmap(i,3))^3);
    end
    
    HH = H'*H; % Assume MF is performed
    
    %% Group_NOMA decision  
        Mhigh   = 0 ;
        Mlow    = 0 ;

        for i = 1 : M        
            if UEmap(i,3) <= d0
                Mhigh = Mhigh + 1 ;             
            else
                Mlow = Mlow + 1 ;
            end 
        end  
        
        Group_H = zeros(1,Mhigh);
        Group_L = zeros(1,Mlow);
        h = 1 ;
        l = 1 ;
        for i = 1 : M        
            if UEmap(i,3) <= d0
                Group_H(h) = i ;
                h = h + 1 ; 
            else
                Group_L(l) = i ;
                l = l + 1 ; 
            end 
        end          
        
        %% Outcome setting 
        P1temp = 0 ;    % H*1 L*1
        P2temp = 0 ;    % H*1 L*0
        P3temp = 0 ;    % H*0 L*1z
        
        %% Two part of RA    
            H_preambleSelect  = zeros(1,Mhigh);
            H_preambleCount   = zeros(1,N);

            L_preambleSelect  = zeros(1,Mlow);
            L_preambleCount   = zeros(1,N);

            for i = 1 : Mhigh
                H_preambleSelect(i) = randi(N);
                H_preambleCount(H_preambleSelect(i)) = H_preambleCount(H_preambleSelect(i)) + 1;
            end

            for i = 1 : Mlow
                L_preambleSelect(i) = randi(N);
                L_preambleCount(L_preambleSelect(i)) = L_preambleCount(L_preambleSelect(i)) + 1;
            end
            Group_H;
            H_preambleSelect;
            H_preambleCount;
            L_preambleSelect;
            L_preambleCount;
            
            Group_NOMA_All = zeros(N,2); 
        
            for i = 1 : N
                if H_preambleCount(i) == 1 && L_preambleCount(i) == 1 
                    P1temp =P1temp + 1; 
                    Group_NOMA_All(i,:) = [Group_H(find(H_preambleSelect == i)) Group_L(find(L_preambleSelect == i))] ;
                elseif H_preambleCount(i) == 1 && L_preambleCount(i) == 0
                    P2temp =P2temp + 1;
                    Group_NOMA_All(i,:) = [Group_H(find(H_preambleSelect == i)) 0] ;
                elseif H_preambleCount(i) == 0 && L_preambleCount(i) == 1
                    P3temp =P3temp + 1;  
                    Group_NOMA_All(i,:) = [0 Group_L(find(L_preambleSelect == i))] ;
                end
            end
    %% 
            Group_NOMA_All;
            c = 1;
            for i = 1 : N
                if Group_NOMA_All(i,1)~=0 ||  Group_NOMA_All(i,2)~=0
                    c = c +1;
                end
            end
            Group_NOMA = zeros(c-1,2);
            d = 1;
            s_1=0;
            s_2=0;
            s_3=0;
            for i = 1:N
                if  Group_NOMA_All(i,1)~=0 ||  Group_NOMA_All(i,2)~=0
                    Group_NOMA(d,:) = Group_NOMA_All(i,:);
                    d=d+1;
                end
                if Group_NOMA_All(i,1)~=0 &&  Group_NOMA_All(i,2)~=0
                    s_1 = s_1 + 1;
                elseif Group_NOMA_All(i,1)~=0 &&  Group_NOMA_All(i,2)==0
                    s_2 = s_2 + 1;
                elseif Group_NOMA_All(i,1)==0 &&  Group_NOMA_All(i,2)~=0
                    s_3 = s_3 + 1;
                end
            end
            
            p_1 =   s_1/N;
            p_2 =   s_2/N;
            p_3 =   s_3/N;
            
            %Group_NOMA
    %%
            P1res(j) = P1temp ;
            P2res(j) = P2temp ;
            P3res(j) = P3temp ;
          
     %% Traditional RACH 
            TR_preambleSelect  = zeros(1,M);
            TR_preambleCount   = zeros(1,N);
            Group_TR = (1:1:M);
            for i = 1 : M
                TR_preambleSelect(i) = randi(N);
                TR_preambleCount(TR_preambleSelect(i)) = TR_preambleCount(TR_preambleSelect(i)) + 1;
            end
    
            GroupTR = zeros(1,N); 
           
            for i = 1 : N
                if TR_preambleCount(i) == 1 
                   GroupTR(i) = Group_TR(find(TR_preambleSelect == i));                
                end
            end       
          
          d = 1;
           for i = 1 : N
               if GroupTR(i)~=0  
                   d = d +1;
               end
           end
           Group_no_zero = zeros(1,d-1);
           e = 1;
           for i = 1:N
                if  GroupTR(i)~=0
                    Group_no_zero(e) = GroupTR(i);
                    e=e+1;
                end
           end    
           Group_OMA= Group_no_zero;        
            
    %% Capacity for NOMA scheme
    for k1 = 1:length(Group_NOMA) % select High UEs
        if (Group_NOMA(k1,1)~=0)
            s_pw = HH(Group_NOMA(k1,1),Group_NOMA(k1,1))^2;
            i_pw = 0;
            if (Group_NOMA(k1,2)~=0)
                i_pw = i_pw + HH(Group_NOMA(k1,1),Group_NOMA(k1,2))'*HH(Group_NOMA(k1,1),Group_NOMA(k1,2));
            end
        else
            s_pw = 0;
            i_pw = 0;
        end
            SINR_H(k1) = s_pw/(i_pw + n_pw);
    end
    
    %
    for k1 = 1:length(Group_NOMA) % select Low UEs
        if (Group_NOMA(k1,2)~=0)
            s_pw = HH(Group_NOMA(k1,2),Group_NOMA(k1,2))^2;
            i_pw = 0;
            if (Group_NOMA(k1,1)~=0)
                i_pw = i_pw + HH(Group_NOMA(k1,1),Group_NOMA(k1,2))'*HH(Group_NOMA(k1,1),Group_NOMA(k1,2));
            end
        else
            s_pw = 0;
            i_pw = 0;
        end
            SINR_L(k1) = s_pw/(i_pw + n_pw);
    end
    
    %
    for k2 = 1:length(Group_NOMA) 
        if (SINR_H(k2) ~= 0 && SINR_L(k2) ~= 0)
            Capacity_H   = log2(1 + Power1*SINR_H(k2)/(1+Power2*SINR_L(k2)));
            Capacity_L   = log2(1 + Power2*SINR_L(k2));
            Capacity(1:2,k2) = [Capacity_H + Capacity_L 1];
        
        elseif (SINR_H(k2) ~= 0 && SINR_L(k2) == 0)
            Capacity(1:2,k2) = [log2(1 + SINR_H(k2)) 2];
        
        elseif (SINR_H(k2) == 0 && SINR_L(k2) ~= 0)
            Capacity(1:2,k2) = [log2(1 + SINR_L(k2)) 3];
        
        end
    end
     %
    for k3 = 1:length(Group_NOMA)
        if (Capacity(2,k3) == 1)
            Capacity_total_temp(k3) = Capacity(1,k3);%*p_1;
        elseif (Capacity(2,k3) == 2)
            Capacity_total_temp(k3) = Capacity(1,k3);%*p_2;
        elseif (Capacity(2,k3) == 3)
            Capacity_total_temp(k3) = Capacity(1,k3);%*p_3;
        end
    end
    sum_Capacity_NOMA(j) = sum(Capacity_total_temp);
    
    %% Capacity of OMA scheme
    count = 1;
    for k1 = Group_OMA % select High UEs
            s_pw = HH(k1,k1)^2;
            i_pw = 0;
            for k2 = Group_OMA
                %if (k2~=k1)
                    %i_pw = i_pw + HH(k1,k2)'*HH(k1,k2);
                %end
            end
            SINR_OMA(count) = s_pw/(i_pw + n_pw);
            count = count + 1;
     end
        sum_Capacity_OMA(j) = sum(log2(1+SINR_OMA));    
    end
    
    
    mean_Capacity_NOMA = mean(sum_Capacity_NOMA);
    mean_Capacity_OMA = mean(sum_Capacity_OMA);

    noma_c(m) = mean_Capacity_NOMA;
    oma_c(m) = mean_Capacity_OMA;
    
    fprintf('loop = %f %%\n', (m*100)/Msize);
end

    
    plot(Mc ,noma_c,'-o');
    hold on;
    plot(Mc ,oma_c,'-*');  
    

