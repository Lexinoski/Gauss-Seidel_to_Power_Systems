clc
clear 
close all
 
fprintf('Aplicação de Gauss Seidel piloto desenvolvida limitada a barras PQ.\n\n')
 
k= input('Fornecer o número de Barras\nk='); %recebe números de barras
%% inicia vetores necessários com valores iniciais em 0;
for i= 1:k;     
    
    E(i)=0;
    Pesp(i)=0;
    Qesp(i)=0;
    Pcalc(i)=0;
    Qcalc(i)=0;
    
    for j=1:k;      %criação das matrizes iniciando 0
        
        Y(i,j)=0;       
        G(i,j)=0;
        B(i,j)=0;
        
    end   
end
 
PV=0;       %contadores dos tipos de barras.
PQ=0;
thetaV=0;
%% recebe todas as informações corretas referente para cada tipo de barra.
for i = 1:k;        
    
    fprintf('\nFornecer o tipo da barra %d: 1-PV ; 2-PQ ; 3-thetaV:',i);
    Tipo(i)= input('');
      
    for j=i:k;      %recebe informações sobre as LTs.
        
        if i~=j;    %matriz auxiliar para armazenamento das impedâncias
            
            fprintf('\nFornecer valor da parte resistiva da LT entre as barras %d e %d em p.u.\nR%d%d=', i, j, i, j);
            R(i,j)= input('');
            fprintf('\nFornecer valor da parte reativa da LT entre as barras %d e %d em p.u.\nX%d%d=', i, j, i, j);
            X(i,j)= input('');
            R(j,i)=R(i,j);
            X(j,i)=X(i,j);
        else
            fprintf('\nFornecer o valor da suceptância capacitiva conectada a barra %d em p.u.\nbsh(%d)=', i, i);
            bsh(i)=input('');
        end
    end
    
    if Tipo(i)==1;
        
        fprintf('\nFornecer o valor de P%d líquida especificada em p.u.\nP%d=',i ,i);
        Pesp(i)=input('');
        fprintf('\nFornecer o valor de V%d em p.u.\nV%d=',i ,i);
        V(i)=input('');  
        PV=PV+1;
   
    elseif Tipo(i)==2;
        
        fprintf('\nFornecer o valor de P%d líquida especificada em p.u.\nP%d=',i ,i);
        Pesp(i)=input('');
        fprintf('\nFornecer o valor de Q%d líquida especificada em p.u.\nQ%d=',i ,i);
        Qesp(i)=input('');
        PQ=PQ+1;
        
    else Tipo(i)==3; 
        
        fprintf('\nFornecer o valor de V%d para barra slack em p.u.\nV%d=',i ,i);
        V(i)=input('');
        fprintf('\nFornecer o valor de theta%d para barra slack em radianos\ntheta%d=',i ,i);
        theta(i)=input('');
        thetaV=thetaV+1;
    end   
end 
Sbase=input('\nAtribua a Potência Base do Sistema em MVA - Sb = ');
 
%% PASSO 1: MONTAR MATRIZ DE ADMITANCIA
for i = 1:k;        %monta matriz de admintância com os dados.
    for j=1:k;
        if i~=j;    %para elementos fora da diagonal principal
        
            Y(i,j)=-(1/(R(i,j)+1i*X(i,j)));
        end
    end
end
 
for i=1:k;                      %para elementos da diagonal principal
    for j=1:k;      
        if i==j;
            for m=1:k;
                if i~=m;        
                    
                    Y(i,i)= Y(i,i)-Y(i,m);    %utiliza elementos ja calculados fora da diagonal que são negativos.
                
                end                
            end
            Y(i,i)= Y(i,i)+1i*bsh(i);  %inclusão do elemento shunt se houver.
        end
    end
end
for i=1:k;              %Abre matriz Y = G+jB;
    for j=1:k;
        
        G(i,j)=real(Y(i,j));
        B(i,j)=imag(Y(i,j));
        
    end
end
%% Identificação de todas a incógnitas e variáveis conhecidas.
fprintf('\nPelo método de Gauss Seidel, serão determinados os seguintes elementos:\n\n');
incog=0;        %contador de incógnitas.
 
for i=1:k;
    if Tipo(i)==1;
        fprintf('theta%d\n',i);
        incog=incog+1;          %acrescenta 1 incógnita.
    elseif Tipo(i)==2;
        fprintf('theta%d\n',i);
        fprintf('V%d\n', i);
        incog=incog+2;          %acrescenta 2 incógnitas.
    end
end  
for i=1:k;
    if Tipo(i)==1;
        fprintf('\nAtribua o kickstart em radianos para theta%d=',i);
        theta(i)=input('');
    elseif Tipo(i)==2;
        fprintf('\nAtribua o kickstart em radianos para theta%d=',i);
        theta(i)=input('');
        fprintf('\nAtribua o kickstart em p.u. para V%d=', i);
        V(i)=input('');
    end
end
Erro=input('\nAtribua o critério de parada Erro<=');
 
for i=1:k;          %Calcula diferença entre Potência específica e calculada.
    if Tipo(i)==1;
        
        DeltaP(i)=Pesp(i)-Pcalc(i);  
   
    elseif Tipo(i)==2;
        
        DeltaP(i)=Pesp(i)-Pcalc(i);
        DeltaQ(i)=Qesp(i)-Qcalc(i);
    end
end
 
%Inicia Vetor de Valores das tensões e potências no modo retangular das iterações;
for i=1:k;
    E(1,i)=V(i)*cos(theta(i))+1i*V(i)*sin(theta(i));
    S(i)=Pesp(i)+1i*Qesp(i);
end
 
%% SOLUÇÃO ITERATIVA.
inter=0; %contador iterações;
DeltaE=1;
while DeltaE>= Erro;
    inter=inter+1;   %Acrescente iteração
    %Valores da tensão e potência modo retangular das iterações;
    for i=1:k;
        if Tipo(i)==3
            E(inter+1,i)=E(inter,i); %Copia valores das barras thetaV
        end
    end
    
    %Operação com a Matriz Y.
    for i=1:k;
        if Tipo(i)==1;
            %NÃO PROGRAMADO
            %NÃO PROGRAMADO
            %NÃO PROGRAMADO
        elseif Tipo(i)==2;
            
            E(inter+1,i)=(1/Y(i,i))*(conj(S(i))/conj(E(inter,i)));
            for j=1:i-1;
                E(inter+1,i)= E(inter+1,i)-(1/Y(i,i))*Y(i,j)* E(inter+1,j);
            end%Developed by Gilberto Lexinoski
            for j=i+1:k;
                E(inter+1,i)= E(inter+1,i)-(1/Y(i,i))*Y(i,j)* E(inter,j);
            end
        end
    end
    for i=1:k;      %aloca valores soluções para barras PQ.
        if Tipo(i)==2;
            V(i)=abs(E(inter+1,i));    %Developed by Gilberto Lexinoski
            theta(i)=atan(imag(E(inter+1,i))/real(E(inter+1,i)));
            %fprintf('\nO novo valor de V%d=',i);
        end
    end
    
    DeltaE=max(abs([E(inter+1,:)]-[E(inter,:)])); %Coleta maior Erro.
end
for i=1:k;          %zera valores de P e Q para recálculo
    for j=1:k;
        Pcalc(i)=0;
        Qcalc(i)=0;
    end
end
    
for i = 1:k;        %Recalcula P e Q de todas as barras.
    for j=1:k;
        
        Pcalc(i)=Pcalc(i)+(V(i)*(V(j)*(G(i,j)*cos(theta(i)-theta(j))+B(i,j)*sin(theta(i)-theta(j)))));
        Qcalc(i)=Qcalc(i)+(V(i)*(V(j)*(G(i,j)*sin(theta(i)-theta(j))-B(i,j)*cos(theta(i)-theta(j)))));
        
    end
end
 fprintf('\nApós %d iterações, se atingiu Erro=%f sendo inferior ao critério de parada E=%f\n', inter, DeltaE, Erro);
fprintf('Se obteve os seguintes valores de convergência:\n');
 
for i=1:k;          %Apresenta os valores de convergência V, theta, P e Q do sistema.
    fprintf('\nV%d = %f p.u.',i,V(i));
    fprintf('\ntheta%d = %f rad',i,theta(i));
    fprintf('\nPcalc%d= %f p.u.    ->    Pcalc%d= %f MW',i,Pcalc(i),i,Pcalc(i)*Sbase);
    fprintf('\nQcalc%d= %f p.u.    ->    Qcalc%d= %f Mvar.\n',i,Qcalc(i),i,Qcalc(i)*Sbase);
end