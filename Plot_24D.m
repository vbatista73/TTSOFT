spec=Plotspec(amp,10,24);
spec(1,:)=1240./spec(1,:);
spec(2,:)=spec(2,:)/max(spec(2,:));
figure;
plot(spec(1,:),spec(2,:));
axis([220 280 0 1.1]);
expt=importdata('ExptYama83.csv');
expt(:,2)=expt(:,2)/max(expt(:,2));
hold on
plot(expt(:,1),expt(:,2));
hold off