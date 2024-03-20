%% plot bubbles' resonance frequencies

colors = {'b','r','k','g','m','c'};
styles = {'-','-.',':','--'};
mrkers = {'s','^','v','o','*','d','x','>','<','+','p','h'};
widths = {0.9,1.1,1.0,1.2};
total_number = 4096;
min_number = 0;
while min_number<total_number
    colors = cat(2,colors,colors);
    styles = cat(2,styles,styles);
    mrkers = cat(2,mrkers,mrkers);
    widths = cat(2,widths,widths);
    min_number = min([length(colors),length(styles),length(mrkers),length(widths)]);
end
colors = colors(1:total_number);
styles = styles(1:total_number);
mrkers = mrkers(1:total_number);
widths = widths(1:total_number);

%% 显示仿真获得的微泡共振频率
model_name = 'Marmottant';
datafile = ['population1_bubble_resonance_',model_name,'.mat'];
load(datafile);

fig = figure; fig.Position = [50,50,1120,480];
subplot(1,2,1);
plt=plot(rads*1e6,Fresn(:,1)/1e6,'b',rads*1e6,Finit(:,1)/1e6,'r'); grid on; hold on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize=4; plt(k).LineWidth=widths{k};
end
title({'Resonance Frequency vs. Bubble Radius',['Model = ', model_name ,', Ambient Overpressure = 0 kPa']});
xlabel('Bubble Radius (\mum)');
ylabel('Frequency (MHz)');
legend('Simulation-derived','Linearization-derived');
subplot(1,2,2);
plt=plot(rads*1e6,Fresn/1e6); grid on; hold on;
for k=1:length(plt)
    plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize=4; plt(k).LineWidth=widths{k};
end
title({'Resonance Frequency vs. Bubble Radius',['Model = ',model_name,', Simulation-derived']});
xlabel('Bubble Radius (\mum)');
ylabel('Frequency (MHz)');
legend(compose("Pov = %2.0f kPa",povs(:)/1e3),'Location','northeast');

