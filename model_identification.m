function [t_min_id,v_min_id, fit] = model_identification(datasets, dIndex, Ts, idx)

% Creo una cartella per conservare le immagini
[ ~, ~ ] = mkdir("Figures/" + datasets(dIndex));

%Caricare un file di condizioni operative
currDataset = load("datasets/" + datasets(dIndex) + ".mat");

%% Sistema da identificare
data = iddata([currDataset.w_mecc_rpm(idx) currDataset.Ce(idx)'], [currDataset.vdRef(idx) currDataset.vqRef(idx) currDataset.Cl(idx)], Ts); 
%data = iddata([ idMis(idx)' iqMis(idx)'], [vdRef(idx) vqRef(idx) theta(idx) Cl(idx)], Ts);
%advice(data)
data = detrend(data);
u = data.u;
y = data.y;
data.OutputName = {'Velocità [RPM]', 'Coppia [N/m]'};
%plot(data)
%data.OutputName = ["w_m_e_c_c_R_P_M" "Cm" ];
%data.InputName = ["vdRef" "vqRef" "Cload"];

%splitting datasets
n = length(currDataset.Cm(idx)');
nt = round(n/2);
train = data(1:nt);
validation = data(nt+1:end);

%% Ordine ottimo 
ord_Max = 30;

% Initialize for better performance
[Jpred_train, Jpred_valid, Jsim_train, Jsim_valid, ...
FPE, AIC, MDL] = deal(zeros(1, ord_Max));


for k = 1:ord_Max
    m2 = arx(train, [k*eye(2) k*ones(2,3) zeros(2,3)]); %modello arx ad ordine k su dati di train
    %m2 = armax(train, [k*eye(2) k*ones(2,3) k*ones(2,1) zeros(2,3)]); %modello ARMAX ad ordine k su dati di train
    %m2 = nlhw(train , [k*eye(2,3) k*ones(2,3) zeros(2,3)]); %modello HW ad ordine k su dati di train
    
    %simulazione del modello 
    y_sim= sim(m2,u); 
    
    %errore di predizione
    ep = pe(m2,data); %prediction error
    err_pred = ep.OutputData;
    
    %errore di simulazione
    err_sim = y-y_sim;
    
    %indici di aderenza (varianza) dell'errore di predizione
    Jpred_train(k) = cov(err_pred(1:nt));
    Jpred_valid(k) = cov(err_pred(nt+1:n));
    
    %indici di aderenza (varianza) dell'errore di simulazione
    Jsim_train(k) = cov(err_sim(1:nt));
    Jsim_valid(k) = cov(err_sim(nt+1:n));
    
    
    %indici di ottimalità
    FPE(k) = fpe(m2); %final prediction error
    AIC(k) = aic(m2);%Akaike information criterion 
    MDL(k) = log(nt)*(k/nt) + log(Jpred_train(k)); %Minimum description length


end

% Find index associated to minimum value of Jpred_valid
[~,v_min_id] = min(Jpred_valid);

% Find index associated to minimum value of Jpred_train
[~,t_min_id] = min(Jpred_train);

... for the following indices, just use the console

%Plotting
x=1:1:ord_Max;

figure (1)
subplot(2,1,1)
plot(x,Jpred_train,  x, Jpred_valid);
xlabel('n - numero di parametri'), ylabel('J predizione'),
legend('J training set', 'J validation set')
title('Andamento degli indici di aderenza in predizione e simulazione')
subplot(2,1,2)
plot(x, Jsim_train, x, Jsim_valid);
xlabel('n - numero di parametri'), ylabel('J simulazione'),
legend('J training set', 'J validation set')
title_string1 = "Indici di aderenza | Dataset: "+ num2str(dIndex);
title(title_string1);
saveas(gcf, get_image_path(datasets(dIndex), "indici_aderenza.png"))


figure (2)
subplot(3,1,1)
plot(x,FPE);
xlabel('n - numero di parametri'), ylabel('FPE'),
title_string2 = "Indici di ottimalità | Dataset: "+ num2str(dIndex);
title(title_string2);
subplot(3,1,2)
plot(x, AIC);
xlabel('n - numero di parametri'), ylabel('AIC'),
subplot(3,1,3)
plot(x, MDL);
xlabel('n - numero di parametri'), ylabel('MDL'),
saveas(gcf, get_image_path(datasets(dIndex), "indice_ottimo.png"))

%% Modello lineare e compare
m = arx(train , [16*eye(2) 16*ones(2,3) zeros(2,3)]);
m1 = arx(train , [20*eye(2) 20*ones(2,3) zeros(2,3)]);

% Other models to try
% m =  armax(train , [16*eye(2) 16*ones(2,3) 1*ones(2,1) zeros(2,3)]);
% m1 = armax(train , [26*eye(2) 26*ones(2,3) 1*ones(2,1) zeros(2,3)]);

% m =  nlhw(train , [15*eye(2,3) 15*ones(2,3) zeros(2,3)]);
% m1 = nlhw(train , [26*eye(2,3) 26*ones(2,3) zeros(2,3)]);

[~,fit_cell,~]=compare(validation, m, m1);
fit = cell2mat(fit_cell);

figure(3)
compare(validation, m, m1)
title_string3 = "Velocità-Coppia | Dataset: "+ num2str(dIndex);
title(title_string3);
saveas(gcf, get_image_path(datasets(dIndex), "model_comparison.png"))

%AutoCorr-XCorr
figure(4)
subplot(211)
resid(m,validation)
subplot(212)
resid(m1,validation)
saveas(gcf, get_image_path(datasets(dIndex), "resid.png"))

clear currDataset data x u y n nt train validation m m1 m2 y_sim ep ...
err_sim err_pred Jpred_train Jpred_valid Jsim_pred Jsim_valid FPE AIC ...
MDL; % Per sicurezza

end

