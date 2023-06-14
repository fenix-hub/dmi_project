function fit_plot(d,vel_fit,torque_fit)
t=1:d;
figure(1)
stem(t,vel_fit(1,:),'b')
hold on 
stem(t,vel_fit(2,:),'--o')
xticks(1:d), xlabel('Datasets'), ylabel('% of fitting'), title('Fitting su Velocità')
legend('m', 'm1')

figure(2)
stem(t,torque_fit(1,:),'b')
hold on 
stem(t,torque_fit(2,:),'--o')
xticks(1:d), xlabel('Datasets'), ylabel('% of fitting'), title('Fitting su Coppia')
legend('m', 'm1')
end
