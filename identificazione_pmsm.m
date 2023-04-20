clc;
clear;

%N.B. Durante l'esecuzione del codice non chiudere le figure. 
Ceidx = 100:50:250;
widx = ["150", "300", "400", "500"];
datasets = strings(0,16);

% Variabili del sistema
Ts = 125;
idx = 5500:15500; % intervallo di regime opportuno
[train_min_id, valid_min_id] = deal(1, zeros(length(datasets)));

i = 1;
for j = 1:length(widx)
    for k = 1:length(Ceidx)
        datasets(i) = "V" + widx(j) + "_" + num2str(Ceidx(k));
        i = i + 1;
    end
end

for d = 1:length(datasets)
    disp("Now working on: " + datasets(d) + "(n°"+num2str(d)+")");
    [t_min_id, v_min_id] = model_identification(datasets, d, Ts, idx);
    valid_min_id(d) = v_min_id;
    train_min_id(d) = t_min_id;
end