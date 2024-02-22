function results = stepwiseRegr(vett,nomiREG,avg_auc,t_onset,DataName1)

% [fid msg]=fopen('REG_TOT2.txt', 'r');
% if(fid>0)
% [vett cont]=fscanf(fid,'%f',[5 450]);
% disp([num2str(cont) ' valori letti...']);
% vett=vett';
% fclose(fid);
% else 
% disp(msg);
% end
% 
% [fid1 msg1]=fopen('regressori2.txt', 'r');
% if(fid1>0)
% nomi5 = textscan(fid1, '%s %s %s %s',1);
% fclose(fid1);
% else 
% disp(msg1);
% end
% 
% 
% for i=1:length(nomi5)
%     nomiREG(i)=cellstr(nomi5{i});
% end
% 
% 
% [fid2 msg2]=fopen('avg_auc2.txt', 'r');
% if(fid2>0)
% [avg_auc cont2]=fscanf(fid2,'%f',[1 450]);
% disp([num2str(cont2) ' valori letti...']);
% avg_auc=avg_auc';
% fclose(fid2);
% else 
% disp(msg2);
% end
% 
% 
% for i=1:length(vett(1,:))
% ind_zeri=find(vett(:,i)==0);
% vett(ind_zeri,:)=[];
% avg_auc(ind_zeri)=[];
% end

% cvr(:,1:(length(vett(1,:))-1))=vett(:,2:end);
cvr=vett;
par=avg_auc;%vett(:,1);


for i=1:length(cvr(1,:))
    H_candidati{i} = cvr(:,i);
end


regressori = nomiREG;
reg= regressori;

n_max=length(H_candidati(1,:));%massimo numero di regressori(colonne di H_candidati)  
n=length(par(:));%numero osservazioni

% H_corrente contiene le colonne che sono presenti in H_canditati(covariate)
H_corrente = ones(n,1); %la prima colonna è di 1 perchè Bo è l'intercetta

beta = (H_corrente'*H_corrente)^-1*H_corrente'*par(:); %trovato ai minimi quadrati (errore normali indipendenti a media nulla e varianza 1 --> sigma_v = Identità)
y_stimato = H_corrente*beta; %calcolato con i parametri appena stimati Y=H0p

% calcolo i residui
residui = par(:) - y_stimato;

SSE = sum(residui.^2);  %variabilita non spiegata
SSY = sum((par(:)-mean(par(:))).^2); %variabilita totale
m = length(H_corrente(1,:))-1;
Ra_corrente = 1-((n-1)/(n-m-1)*SSE/SSY); %% score corrente
%Ra_corr=1-SSE/SSY;

%il primo regressore sarà sempre l'intercetta
regressori_scelti = ['intercept'];

%%Scelta del modello di regressione mediante procedura Stepwise forward
 
for j = 1 : n_max %deve ciclare su tutti i regressori che considero
   
    %disp(['Ra corrente = ',num2str(Ra_corr)])
    
    for i = 1:length(H_candidati(1,:))
        %aggiungo di volta in volta un regressore e controllo lo score del
        %modello--> compromesso tra numero di iterazioni dovute a
        %combinazioni numerosissime e affidabilità modello (li prendo "ad uno ad uno")
        
        H = [H_corrente cell2mat(H_candidati(:,i))];
        beta = (H'*H)^-1*H'*par(:);
        y_stimato = H*beta; % Y calcolato con i parametri appena stimati
    
        % calcolo i residui
        residui = par(:) - y_stimato; 
        
        SSE = sum(residui.^2); %(Yoss-Ypred)^2 variabilità non spiegata
        SSY = sum((par(:)-mean(par(:))).^2);%(Ypred-Ymed)^2 variabilità totale
        
        m = length(H(1,:))-1;%#regressori correnti
        Ra(i) = 1-((n-1)/(n-m-1)*SSE/SSY); %% score i-esimo (R aggiustato)
         
        %R(i) = 1-SSE/SSY;
     end
    
   %Ra
    
   [r_max imax] = max(Ra); %cerco lo score più alto raggiunto nel set di regressori
    if(r_max>Ra_corrente)%se il punteggio si è alzato rispetto a quello corrente allora lo aggiungo, altrimenti no
        % H_corrente contenva solo una colonna di 1, ora aggiungo il migliore trovato sopra, quindi riparte il ciclo for su i con
        % H_corrente costituita da due colonne e in H_candidati cancello la colonna del massimo regressore trovato
        H_corrente = [H_corrente cell2mat(H_candidati(:,imax))];
        H_candidati(:,imax) = [];
        Ra_corrente = Ra(imax);
        index{i}=imax;
        
        regressori_scelti = [regressori_scelti regressori(imax)];    
       
        regressori(imax) = [];%annullo la variabile per renderla accessibile alla prossima iterazione
        
    else
        break %quando lo score si "abbassa" esco dal ciclo
    end
end
Ra=[];
%R_corr=[];
beta= (H_corrente'*H_corrente)^-1*H_corrente'*par(:);
SSE = sum(residui.^2);
SSY = sum((par(:)-mean(par(:))).^2);
disp(['Selected regressors: ', regressori_scelti])
%solo per farmeli stampare
Beta=beta
%per ciascuno dei miei parametri Vd, t_short, t_long, FRA mi assoccia i beta migliori trovati
B=beta;

%%Test d'ipotesi con H0 = B0,B1,B2,...=0 livello di fiducia al 99%
%fatto nel ciclo z
SSR =  SSY - SSE;
F = (SSR/length(beta)) / (SSE/(n-length(beta)-1))
% FINV   Inverse of the F cumulative distribution function.
%     X=FINV(P,V1,V2) returns the inverse of the F distribution 
%     function with V1 and V2 degrees of freedom, at the values in P.
val_cr_F = finv(0.99,length(beta),n-length(beta)-1)
if F > val_cr_F
    disp('F test result: rejection H0 = beta(k-th) is different from zero')
else
    disp('F test result: I cannot reject H0 =beta(k-th) is equal to zero')
end 


%%Costruiamo gli intervalli di confidenza dei parametri usando la t-student con livello di fiducia al 99% -->
%%alfa/2=0.005
% si può dimostrare che il valor medio mu_beta di ognuno degli stimatori
% ottenuti con il metodo dei minimi quadrati beta_stimato è proprio il
% parametro beta(i)
%la varianza s^2 di beta_stimato si ottiene moltiplicando s2 per una
%funzione data dall'elemento sulla diagonale principale della matrice
%(H_corrente'*H_corrente)^-1

var_beta = inv(H_corrente'*H_corrente);%varianza dei beta
s2 = SSE / ( n-length(beta) );%varianza delle stime

for i = 1 : length(beta)
    
    s_beta_i(i) = sqrt( var_beta(i,i)*s2);
        %estremi superiore ed inferiore
    b_inf(i) = beta(i) - ( s_beta_i(i) * tcdf(0.005,n-length(beta)));
    b_sup(i) = beta(i) + ( s_beta_i(i) * tcdf(0.005,n-length(beta)));
    
end
   Binf=b_inf;
   Bsup=b_sup;
   
   b_inf=[];
   b_sup=[];

%calcolo il coefficiente di variazione facendo il rapporto tra deviazione standard e valori
matrice_cov = diag( (H_corrente'*H_corrente)^-1 );
disp('Coefficients of variation of parameters:')
CV = sqrt(matrice_cov)./Beta


disp(' ')
disp('Parameters estimated with my data:')

results.y_estimates = y_stimato;
results.par = par;
results.reg = reg;
results.sel_reg = regressori_scelti;
results.beta = Beta;

% myfilenameWRITE= sprintf([num2str(t_onset) '_' DataName1 '.xlsx'])
% 
%    xlswrite(myfilenameWRITE,y_stimato,'Sheet1','c2');
%    xlswrite(myfilenameWRITE,par,'Sheet1','b2');
%    xlswrite(myfilenameWRITE,reg','Sheet1','e2')
%    xlswrite(myfilenameWRITE,regressori_scelti','Sheet1','g2')
%    xlswrite(myfilenameWRITE,Beta,'Sheet1','h2')
   
   

  %BlandAltman(par, y_stimato)
  
  %sgtitle([' Y-predicted: ' '____________' DataName1 ' ' ' ' ' ' ' Regressors: ' '____________' regressori_scelti])

  
   
% for i=1:length(par)
%     if par(i)>=0.11;
%         malann{i}='HIGH_MSNA';
%     else
%         malann{i}='LOW_MSNA';
%     end
% end
% 
% 
% 
% [X,Y,T,AUC] = perfcurve(malann,y_stimato,'HIGH_MSNA'); 
% 
% AUC
% 
% % Plot the ROC curve.
% figure
% 
% for id_t=1:length(y_stimato)
% distance(id_t)= sqrt((X(id_t)^2+(Y(id_t)-1)^2));
% end
% 
% [sss tmin]=min(distance);
%         
% fill_color = [11/255, 208/255, 217/255];
% fill([X; 1], [Y; 0], fill_color,'FaceAlpha',0.5);
% hold on; plot(X, Y, '-b', 'LineWidth', 2);
% hold on; plot(X(tmin), Y(tmin), 'or', 'MarkerSize', 10);
% hold on; plot(X(tmin), Y(tmin), 'xr', 'MarkerSize', 12);
% hold off; axis square; grid on; grid minor; xlabel('1 - specificity'); ylabel('sensibility');
% title('ROC for Classification by Logistic Regression')
% legend(['AUC = ' num2str(AUC)]);
% 
% 
% a=sprintf('T = %.3f -->  ',T(tmin));
% text(X(tmin), Y(tmin),a,'HorizontalAlignment','right')


end