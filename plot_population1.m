% close all;

%% 画图设置项
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

%% 选择项设置
POV_FOC = [0,10,25] * 1e3;
[~,POV_LST] = min(abs(POV_FOC - povs(:)));
FRQ_FOC = [2.0,2.5,3.0,3.5,4.0] * 1e6;
[~,FRQ_LST] = min(abs(FRQ_FOC - frqs(:)));
PAC_FOC = [100,300,500] * 1e3;
[~,PAC_LST] = min(abs(PAC_FOC - pacs(:)));

%% 数据处理
xsh = zeros(length(rads),length(povs),length(frqs),length(pacs));
xfd = xsh; xuh = xsh; xhm = xsh;
MAX_EXPAND = xsh; MAX_COMPRES = xsh;

for pov_ndx = 1:length(povs)
    Pov = povs(pov_ndx);
    for frq_ndx = 1:length(frqs)
        Frq = frqs(frq_ndx);
        for pac_ndx = 1:length(pacs)
            Pac = pacs(pac_ndx);
            
            Tt = Tnodes{pov_ndx,frq_ndx,pac_ndx};
            Pt = Pnodes{pov_ndx,frq_ndx,pac_ndx};
            Ti = Tincident{pov_ndx,frq_ndx,pac_ndx};
            Ri = Rincident{pov_ndx,frq_ndx,pac_ndx};
            
            for rad_ndx = 1:length(rads)
                Rad = rads(rad_ndx);
                
                %% 信号处理流程
                % 提取振动信号
                timein = Ti;
                radius = Ri(:,rad_ndx);
                TIN_LIST = timein<=((sg_window_end-sg_window_beg)+timein(1));
                
                % 提取回声信号
                timesg = Tt(:,rad_ndx);
                signal = Pt(:,rad_ndx);
                ROI_LIST = (timesg>=sg_window_beg) & (timesg<=sg_window_end);
                
                % 振动信号分析
                MAX_EXPAND(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = max(radius) / Rinit(rad_ndx);
                MAX_COMPRS(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = min(radius) / Rinit(rad_ndx);
                
                % ROI信号的分析
                roi_signal = signal(ROI_LIST);
                [pows,freq,powsc] = periodogram(roi_signal,...
                    rectwin(length(roi_signal)),length(roi_signal),probe_fs,...
                    'onesided','power','ConfidenceLevel',0.95);
                resp = 10*log10(pows*2);
                respc = 10*log10(powsc*2);
                
                % 谐波幅度的计算
                fr_wd = 0.25e6; %Hz
%                 xsh(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = mean(resp(freq<Frq*0.5+fr_wd & freq>=Frq*0.5-fr_wd));
%                 xfd(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = mean(resp(freq<Frq*1.0+fr_wd & freq>=Frq*1.0-fr_wd));
%                 xuh(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = mean(resp(freq<Frq*1.5+fr_wd & freq>=Frq*1.5-fr_wd));
%                 xhm(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = mean(resp(freq<Frq*2.0+fr_wd & freq>=Frq*2.0-fr_wd));
                xsh(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = interp1(freq(freq<Frq*0.5+fr_wd & freq>=Frq*0.5-fr_wd),resp(freq<Frq*0.5+fr_wd & freq>=Frq*0.5-fr_wd),Frq*0.5);
                xfd(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = interp1(freq(freq<Frq*1.0+fr_wd & freq>=Frq*1.0-fr_wd),resp(freq<Frq*1.0+fr_wd & freq>=Frq*1.0-fr_wd),Frq*1.0);
                xuh(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = interp1(freq(freq<Frq*1.5+fr_wd & freq>=Frq*1.5-fr_wd),resp(freq<Frq*1.5+fr_wd & freq>=Frq*1.5-fr_wd),Frq*1.5);
                xhm(rad_ndx,pov_ndx,frq_ndx,pac_ndx) = interp1(freq(freq<Frq*2.0+fr_wd & freq>=Frq*2.0-fr_wd),resp(freq<Frq*2.0+fr_wd & freq>=Frq*2.0-fr_wd),Frq*2.0);

                
                %% 显示不同大小微泡的振动、回声及其频谱
                if length(rads)>=1 && length(povs)==1 && length(frqs)==1 && length(pacs)==1
                    
                    fig10000=figure(10000); fig10000.Position=[880 140 960 960];
                    
                    subplot(length(rads),3,3*(rad_ndx-1)+1);
                    plot(timein*1e6,radius/Rinit(rad_ndx,pov_ndx),'LineWidth',1);
                    xlabel('Time (\mus)'); ylabel('R/R0');
                    title_str1 = ['\itR_{0} = ',num2str(Rinit(rad_ndx,pov_ndx)*1e6,'%2.2f'),...
                        ' \mum, f_{r} = ',num2str(Finit(rad_ndx,pov_ndx)/1e6,'%2.2f MHz')];
                    title_str2 = ['\itE = ',num2str(max(radius)/Rinit(rad_ndx,pov_ndx),'%2.2f'),...
                        ', C = ',num2str(min(radius)/Rinit(rad_ndx,pov_ndx),'%2.2f'),...
                        ', E/C = ',num2str((max(radius)-Rinit(rad_ndx,pov_ndx))/(Rinit(rad_ndx,pov_ndx)-min(radius)),'%1.2f')];
                    if rad_ndx==1
                        title({'Radius Dynamics','',title_str1,title_str2});
                    else
                        title({title_str1,title_str2});
                    end
                    grid on; hold on;
                    YL = [0.5,2] ;
                    if strcmp(bubble_model,'Marmottant')
                        YL = [0.5,2];
                    end
                    ti_roi = timein(TIN_LIST);
                    rectangle('Position',[ti_roi(1)*1e6,YL(1),ti_roi(end)*1e6-ti_roi(1)*1e6,YL(2)-YL(1)],'EdgeColor','g','LineStyle','-','LineWidth',1)
                    ylim(YL); xlim([ti_roi(1)*1e6-1,ti_roi(end)*1e6+3]);
                    
                    h=subplot(length(rads),3,3*(rad_ndx-1)+2);
                    plot(timesg*1e6,signal,'LineWidth',1);  hold on;
                    xlabel('Time (\mus)'); ylabel('Pressure (Pa)');
                    title_str = ['\itPNP = ',num2str(min(signal),'%2.2f'),' Pa'];
                    if rad_ndx==1
                        title({'Scattered Acoustic Signal','','',title_str});
                    else
                        title({title_str});
                    end
                    grid on; hold on;
                    YL = [-80,100] ;
                    if strcmp(bubble_model,'Marmottant')
                        YL = [-50,50];
                    end
                    tr_roi = timesg(ROI_LIST);
                    rectangle('Position',[tr_roi(1)*1e6,YL(1),tr_roi(end)*1e6-tr_roi(1)*1e6,YL(2)-YL(1)],'EdgeColor','g','LineStyle','-','LineWidth',1)
                    ylim(YL); xlim([tr_roi(1)*1e6-1,tr_roi(end)*1e6+1]);
                    
                    h=subplot(length(rads),3,3*(rad_ndx-1)+3);
                    plot(freq/1e6,resp, 'LineWidth',1);
                    grid on; hold on;
                    xlim([0,10]);
                    h.XTick = 0:(Frq/1e6/2):10;
                    xlabel('Frequency (MHz)'); ylabel('Magnitude (dB)');
                    title_str1 = ['\itSH = ',num2str(xsh(rad_ndx,pov_ndx,frq_ndx,pac_ndx),'%2.1f'),...
                        ', H1 = ',num2str(xfd(rad_ndx,pov_ndx,frq_ndx,pac_ndx),'%2.1f'),...
                        ', H2 = ',num2str(xhm(rad_ndx,pov_ndx,frq_ndx,pac_ndx),'%2.1f dB')];
                    if rad_ndx==1
                        title({'Magnitude Spectrum of Acoustic Signal','','',title_str1});
                    else
                        title(title_str1);
                    end
                    YL = [-40,30];
                    if strcmp(bubble_model,'Marmottant')
                        YL = [-40,30];
                    end
                    ylim(YL);
                end
                
                %% 显示在不同声压作用下的不同大小微泡的振动、回声及其频谱
                if length(rads)>1 && length(povs)==1 && length(frqs)==1 && length(pacs)>1
                    lgds_str = compose('Pac = %2.0f kPa',pacs(PAC_LST)/1e3);
                    if ismember(pac_ndx,PAC_LST)
                        fig10030=figure(10030); fig10030.Position = [100,100,1680,960];
                        h=subplot(length(rads),3,(rad_ndx-1)*3+1);
                        plt=plot(timein*1e6,radius/Rinit(rad_ndx,pov_ndx),'LineWidth',1); grid on; hold on;
                        k = pac_ndx;
                        plt.Color=colors{k}; plt.LineStyle=styles{k}; plt.Marker='none'; plt.LineWidth=widths{k};
                        if pac_ndx==length(pacs)
                            legend(lgds_str,'Location','best');
                            xlabel('Time (\mus)'); ylabel('R/R0');
                            title({['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.2f'),' \mum'],...
                                ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'Ambient Overpressure = ',num2str(Pov/1e3),' kPa'],...
                                ['Radius Dynamics']}); grid on; hold on; hold on;
                            title({'Radius Dynamics',['\itR_{0} = ', num2str(R00(rad_ndx)*1e6,'%2.2f'),' \mum, f_{r} = ', num2str(Finit(rad_ndx,pov_ndx)/1e6,'%2.2f'),' MHz']});
                            YL = [0.4,2.8];
                            if strcmp(bubble_model,'Marmottant')
                                YL = [0.4,2];
                            end
                            ylim(YL);
                            ti_roi = timein(TIN_LIST);
                            xlim([ti_roi(1)*1e6-1,ti_roi(end)*1e6+1]);
                            rectangle('Position',[ti_roi(1)*1e6,YL(1),ti_roi(end)*1e6-ti_roi(1)*1e6,YL(2)-YL(1)],'EdgeColor','g','LineStyle','-','LineWidth',1)
                        end
                        
                        h=subplot(length(rads),3,(rad_ndx-1)*3+2);
                        plt=plot(timesg*1e6,signal,'LineWidth',1); grid on; hold on;
                        k = pac_ndx;
                        plt.Color=colors{k}; plt.LineStyle=styles{k}; plt.Marker='none'; plt.LineWidth=widths{k};
                        if pac_ndx==length(pacs)
                            legend(lgds_str,'Location','best');
                            xlabel('Time (\mus)'); ylabel('Pressure (Pa)');
                            title({'Scattered Acoustic Signal',''}); grid on; hold on; hold on;
                            YL = [-100,160];
                            if strcmp(bubble_model,'Marmottant')
                                YL = [-100,100];
                            end
                            ylim(YL);
                            tr_roi = timesg(ROI_LIST);
                            xlim([tr_roi(1)*1e6-1,tr_roi(end)*1e6+1]);
                            rectangle('Position',[tr_roi(1)*1e6,YL(1),tr_roi(end)*1e6-tr_roi(1)*1e6,YL(2)-YL(1)],'EdgeColor','g','LineStyle','-','LineWidth',1)
                        end
                        
                        h=subplot(length(rads),3,(rad_ndx-1)*3+3);
                        plt=plot(freq/1e6,resp-20*log10(Pac), 'LineWidth',1); grid on; hold on;
                        k = pac_ndx;
                        plt.Color=colors{k}; plt.LineStyle=styles{k}; plt.Marker='none'; plt.LineWidth=widths{k};
                        if pac_ndx==length(pacs)
                            xlabel('Frequency (MHz)'); ylabel('Amplitude (dB)');
                            title({'Relative Spectrum of Scattered Acoustic Signal',''});
                            h.XTick = 0:(Frq/1e6/2):10;
                            xlim([Frq/1e6/2-0.5, Frq/1e6+1]);
                            xlim([0,7.5])
                            legend(lgds_str,'Location','best');
                        end
                    end
                end
                
                %% 显示在不同环境压力下的不同大小微泡的振动、回声及其频谱
                if length(rads)>1 && length(rads)<=7 && length(povs)>1 && length(frqs)==1 && length(pacs)==1
                    lgds_str = compose('Pov = %2.0f kPa',povs/1e3);
                    fig10040=figure(10040); fig10040.Position = [100,100,1680,960];
                    h=subplot(length(rads),3,(rad_ndx-1)*3+1);
                    plt=plot(timein*1e6,radius/Rinit(rad_ndx,pov_ndx),'LineWidth',1); grid on; hold on;
                    k = pov_ndx;
                    plt.Color=colors{k}; plt.LineStyle=styles{k}; plt.Marker='none'; plt.LineWidth=widths{k};
                    if pov_ndx==length(povs)
                        legend(lgds_str,'Location','best');
                        xlabel('Time (\mus)'); ylabel('R/R0');
                        if pulse_negative==1
                            title({['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.2f'),' \mum'],...
                                ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'Pulse Magnitude = -',num2str(Pac/1e3),' kPa'],...
                                ['Radius Dynamics']}); grid on; hold on; hold on;
                        else
                            title({['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.2f'),' \mum'],...
                                ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'Pulse Magnitude = ',num2str(Pac/1e3),' kPa'],...
                                ['Radius Dynamics']}); grid on; hold on; hold on;
                        end
                        title({'Radius Dynamics',['\itR_{0} = ', num2str(R00(rad_ndx)*1e6,'%2.2f'),' \mum, f_{r} = ', num2str(Finit(rad_ndx,1)/1e6,'%2.2f'),' MHz']});
                        YL = [0.5,2.0];
                        if strcmp(bubble_model,'Marmottant')
                            YL = [0.5,2.0];
                        end
                        ylim(YL);
                        ti_roi = timein(TIN_LIST);
                        xlim([ti_roi(1)*1e6-1,ti_roi(end)*1e6+1]);
                        rectangle('Position',[ti_roi(1)*1e6,YL(1),ti_roi(end)*1e6-ti_roi(1)*1e6,YL(2)-YL(1)],'EdgeColor','g','LineStyle','-','LineWidth',1)
                    end
                    
                    h=subplot(length(rads),3,(rad_ndx-1)*3+2);
                    plt=plot(timesg*1e6,signal,'LineWidth',1); grid on; hold on;
                    k = pov_ndx;
                    plt.Color=colors{k}; plt.LineStyle=styles{k}; plt.Marker='none'; plt.LineWidth=widths{k};
                    if pov_ndx==length(povs)
                        legend(lgds_str,'Location','best');
                        xlabel('Time (\mus)'); ylabel('Pressure (Pa)');
                        title({'Received Acoustic Signal',''}); grid on; hold on; hold on;
                        YL = [-80,80];
                        if strcmp(bubble_model,'Marmottant')
                            YL = [-60,60];
                        end
                        ylim(YL);
                        tr_roi = timesg(ROI_LIST);
                        xlim([tr_roi(1)*1e6-1,tr_roi(end)*1e6+1]);
                        rectangle('Position',[tr_roi(1)*1e6,YL(1),tr_roi(end)*1e6-tr_roi(1)*1e6,YL(2)-YL(1)],'EdgeColor','g','LineStyle','-','LineWidth',1)
                    end
                    
                    h=subplot(length(rads),3,(rad_ndx-1)*3+3);
                    plt=plot(freq/1e6,resp, 'LineWidth',1); grid on; hold on;
                    k = pov_ndx;
                    plt.Color=colors{k}; plt.LineStyle=styles{k}; plt.Marker='none'; plt.LineWidth=widths{k};
                    if pov_ndx==length(povs)
                        xlabel('Frequency (MHz)'); ylabel('Amplitude (dB)');
                        title({'Spectrum of Received Acoustic Signal',''});
                        h.XTick = 0:(Frq/1e6/2):10;
                        xlim([Frq/1e6/2-0.5, Frq/1e6+1]);
                        xlim([0,7.5]);
                        YL = [-40,30];
                        if strcmp(bubble_model,'Marmottant')
                            YL = [-40,30];
                        end
                        ylim(YL);
                        legend(lgds_str,'Location','best');
                    end
                    
                end
                
            end
        end
    end
end

%% 显示不同大小微泡的最大膨胀和收缩，以及次谐波幅值
if length(rads)>1 && length(povs)==1 && length(frqs)==1 && length(pacs)==1
    
    fig=figure(10010); fig.Position = [680 40 560 480];
    plt=plot(rads(:)*1e6,[xsh(:,pov_ndx,frq_ndx,pac_ndx),...
        xfd(:,pov_ndx,frq_ndx,pac_ndx),...
        xuh(:,pov_ndx,frq_ndx,pac_ndx),...
        xhm(:,pov_ndx,frq_ndx,pac_ndx)],'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend('Subharmonic','Fundamental','Ultraharmonic','2^{nd} Harmonic','Location','best');
    xlabel('Bubble Radius(\mum)'); ylabel('Amplitude(dB)');
    if pulse_negative==1
        title({'Bubble Harmonics vs. Bubble Radius',...
            ['Model = ', bubble_model, ', Pulse Magnitude = -', num2str(Pac/1e3,'%2.0f'),' kPa'],...
            ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'OverPressure = ',num2str(Pov/1e3),' kPa']});
    else
        title({'Bubble Harmonics vs. Bubble Radius',...
            ['Model = ', bubble_model, ', Pulse Magnitude = ', num2str(Pac/1e3,'%2.0f'),' kPa'],...
            ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'OverPressure = ',num2str(Pov/1e3),' kPa']});
    end
    xlim([rads(1),rads(end)]*1e6);
    
    fig=figure(10020); fig.Position = [180 40 1120 480];
    subplot(1,2,1);
    plt=plot(rads*1e6,[MAX_EXPAND(:),MAX_COMPRS(:)],'-*','LineWidth',1); grid on; hold on;
    xlabel('Bubble Radius (\mum)'); ylabel('E&C (100%)');
    if pulse_negative==1
        title({'Maximum Expansion & Compression',...
            ['Model = ', bubble_model, ', Pulse Magnitude = -', num2str(Pac/1e3,'%2.2f'),' kPa'],...
            ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'OverPressure = ',num2str(Pov/1e3),' kPa']});
    else
        title({'Maximum Expansion & Compression',...
            ['Model = ', bubble_model, ', Pulse Magnitude = ', num2str(Pac/1e3,'%2.2f'),' kPa'],...
            ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'OverPressure = ',num2str(Pov/1e3),' kPa']});
    end
    legend('Maximum Expansion','Maximum Compression','Location','best');
    xlim([rads(1),rads(end)]*1e6);
    subplot(1,2,2);
    plt=plot(rads*1e6,xsh(:,pov_ndx,frq_ndx,pac_ndx)); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend('Subharmonic','Location','best');
    xlabel('Bubble Radius(\mum)'); ylabel('Amplitude(dB)');
    if pulse_negative==1
        title({'Bubble Subharmonics vs. Bubble Radius',...
            ['Model = ', bubble_model, ', Pulse Magnitude = -', num2str(Pac/1e3,'%2.0f'),' kPa'],...
            ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'OverPressure = ',num2str(Pov/1e3),' kPa']});
    else
        title({'Bubble Subharmonics vs. Bubble Radius',...
            ['Model = ', bubble_model, ', Pulse Magnitude = ', num2str(Pac/1e3,'%2.0f'),' kPa'],...
            ['Pulse Frequency = ',num2str(Frq/1e6),' MHz, ', 'OverPressure = ',num2str(Pov/1e3),' kPa']});
    end
    xlim([rads(1),rads(end)]*1e6);
%         figure(10020);
%         subplot(1,2,1); legend('Negative Pulse - Maximum Expansion','Negative Pulse - Maximum Compression','Positive Pulse - Maximum Expansion','Positive Pulse - Maximum Compression','Location','best');
%         subplot(1,2,2); legend('Negative Pulse','Positive Pulse','Location','best');
end

%% 显示不同大小微泡的次谐波幅度-环境压力的关系
if length(rads)>=1 && length(povs)>1 && length(frqs)==1 && length(pacs)==1
    fig10050=figure(10050); fig10050.Position = [50 500 560 480];
    X = povs(:) / 1e3; Y = xsh(:,:,frq_ndx,pac_ndx)';
    plt=plot(X,Y,'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Overpressure (kPa)'); ylabel('Subharmonic Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = -',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    else
        title({'Subharmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = ',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    end
    
end

%% 显示不同大小微泡的次谐波幅度-脉冲强度的关系
if length(rads)>1 && length(povs)==1 && length(frqs)==1 && length(pacs)>1
    fig10090=figure(10090); fig10090.Position = [80 540 560 480];
    X = pacs(:) / 1e3; Y = xsh(:,pov_ndx,frq_ndx,:); Y = reshape(Y,length(rads),length(pacs)); Y = Y';
    plt=plot(X,Y,'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).MarkerSize=4; plt(k).LineWidth=1.0;
    end
    legend(compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Pulse Magnitude (kPa)'); ylabel('Subharmonic Amplitude (dB)');
    title({'Subharmonic Amplitude vs. Pulse Magnitude',...
        ['Model = ', bubble_model],...
        ['Ambient Overpressure = ',num2str(Pov/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
end

%% 显示不同大小微泡的次谐波幅度-发射频率的关系
if length(rads)>1 && length(povs)==1 && length(frqs)>1 && length(pacs)==1
    
    fig10100=figure(10100); fig10100.Position = [550 50 560 480];
    X = frqs(:)./(ones(length(frqs),1)*Finit(:,pov_ndx)'); Y = reshape(xsh(:,1,:),length(rads),length(frqs)); Y = Y';
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize = 3; plt(k).LineWidth=widths{k};
    end
    legend(compose("R0 = %2.1f \x03BCm",rads(:)*1e6),'Location','best');
    xlabel('Standardized Frequency (100%)'); ylabel('Subharmonic Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Pulse Frequency',...
            ['Model = ', bubble_model, ', Ambient Overpressure = ', num2str(povs(pov_ndx)/1e3,'%2.1f'),' kPa'],...
            ['Pulse Magnitude = -',num2str(pacs(pac_ndx)/1e3),' kPa']});
    else
        title({'Subharmonic Amplitude vs. Pulse Frequency',...
            ['Model = ', bubble_model, ', Ambient Overpressure = ', num2str(povs(pov_ndx)/1e3,'%2.1f'),' kPa'],...
            ['Pulse Magnitude = ',num2str(pacs(pac_ndx)/1e3),' kPa']});
    end
    xlim([0,4]);
    
end

%% 显示不同环境压力下的次谐波幅度-发射频率的关系
if length(rads)==1 && length(povs)>1 && length(frqs)>1 && length(pacs)==1
    fig10060=figure(10060); fig10060.Position = [50 50 560 480];
    X = frqs(:)/Finit(rad_ndx,1); Y = reshape(xsh(1,:,:),length(povs),length(frqs)); Y = Y(POV_LST,:)';
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize = 3; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.0f kPa",povs(POV_LST)/1e3),'Location','northeast');
    xlabel('Standardized Frequency (f / f_{r})'); ylabel('Subharmonic Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Standardized Frequency',...
            ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
            ['Pulse Magnitude = -',num2str(pacs(pac_ndx)/1e3),' kPa']});
    else
        title({'Subharmonic Amplitude vs. Standardized Frequency',...
            ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
            ['Pulse Magnitude = ',num2str(pacs(pac_ndx)/1e3),' kPa']});
    end
    xlim([frqs(1)/Finit(rad_ndx,1),frqs(end)/Finit(rad_ndx,1)]);
    
    fig10070=figure(10070); fig10070.Position = [550 50 560 480];
    X = frqs(:)/1e6; Y = reshape(xsh(1,:,:),length(povs),length(frqs)); Y = Y(POV_LST,:)';
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize = 3; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.0f kPa",povs(POV_LST)/1e3),'Location','northeast');
    xlabel('Pulse Frequency (MHz)'); ylabel('Subharmonic Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Standardized Frequency',...
            ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
            ['Pulse Magnitude = -',num2str(pacs(pac_ndx)/1e3),' kPa']});
    else
        title({'Subharmonic Amplitude vs. Standardized Frequency',...
            ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
            ['Pulse Magnitude = ',num2str(pacs(pac_ndx)/1e3),' kPa']});
    end
    xlim([frqs(1)/1e6,frqs(end)/1e6]);
    
end

%% 显示不同环境压力下的次谐波幅度-脉冲强度的关系
if length(rads)==1 && length(povs)>1 && length(frqs)==1 && length(pacs)>1
    fig10110=figure(10110); fig10110.Position = [50 50 560 480];
    X = pacs(:)/1e3; Y = reshape(xsh(1,:,1,:),length(povs),length(pacs)); Y = Y(POV_LST,:)';
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize = 3; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.0f kPa",povs(POV_LST)/1e3),'Location','northeast');
    xlabel('Pulse Magnitude (kPa)'); ylabel('Subharmonic Amplitude (dB)');
    title({'Subharmonic Amplitude vs. Pulse Magnitude',...
        ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
        ['Pulse Frequency = ',num2str(frqs(frq_ndx)/1e6),' MHz']});
    xlim([pacs(1)/1e3,pacs(end)/1e3]);
    
    fig10111=figure(10111); fig10111.Position = [50 50 560 480];
    
    X = pacs(:)/1e3; Y = reshape(xsh(1,:,1,:),length(povs),length(pacs)) - ones(length(povs),1) * 20 * log10(pacs); Y = Y(POV_LST,:)';
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize = 3; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.0f kPa",povs(POV_LST)/1e3),'Location','northeast');
    xlabel('Pulse Magnitude (kPa)'); ylabel('Relative Subharmonic Amplitude (dB)');
    title({'Relative Subharmonic Amplitude vs. Pulse Magnitude',...
        ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
        ['Pulse Frequency = ',num2str(frqs(frq_ndx)/1e6),' MHz']});
    xlim([pacs(1)/1e3,pacs(end)/1e3]);
end

%% 显示不同环境压力下的次谐波幅度-微泡大小的关系
if length(rads)>1 && length(povs)>1 && length(frqs)==1 && length(pacs)==1
    fig10120=figure(10120); fig10120.Position = [50 50 560 480];
    X = rads(:)*1e6; Y = reshape(xsh(:,:,1,1),length(rads),length(povs)); Y = Y(:,POV_LST);
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker='o'; plt(k).MarkerSize = 3; plt(k).LineWidth=widths{k};
    end
    legend(compose("Pov = %2.0f kPa",povs(POV_LST)/1e3),'Location','northeast');
    xlabel('Bubble Radius (\mum)'); ylabel('Subharmonic Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Bubble Radius',...
            ['Model = ', bubble_model, ', Pulse Magnitude = -', num2str(pacs(pac_ndx)/1e3,'%2.0f'),' kPa'],...
            ['Pulse Frequency = ',num2str(frqs(frq_ndx)/1e6),' MHz']});
    else
        title({'Subharmonic Amplitude vs. Bubble Radius',...
            ['Model = ', bubble_model, ', Pulse Magnitude = ', num2str(pacs(pac_ndx)/1e3,'%2.0f'),' kPa'],...
            ['Pulse Frequency = ',num2str(frqs(frq_ndx)/1e6),' MHz']});
    end
    xlim([rads(1)*1e6,rads(end)*1e6]);
    
end

%% 显示不同发射频率的次谐波幅度-环境压力的关系
if length(rads)==1 && length(povs)>1 && length(frqs)>1 && length(pacs)==1
    fig10080=figure(10080); fig10080.Position = [50 550 560 480];
    X = povs(:)/1e3; Y = reshape(xsh(1,:,:),length(povs),length(frqs)); Y = Y(:,FRQ_LST);
    plt = plot(X,Y,'Marker','o'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    %compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best'
    legend(compose("f = %2.1f MHz, \x0394SH = %2.2f dB",frqs(FRQ_LST)'/1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Ambient Overpressure (kPa)'); ylabel('Subharmonic Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Standardized Frequency',...
            ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
            ['Pulse Magnitude = -',num2str(pacs(pac_ndx)/1e3),' kPa'],...
            ['Reference Rosonance Frequency f_{r} = ',num2str(Finit(rad_ndx,1)/1e6,'%2.2f'),' MHz']});
    else
        title({'Subharmonic Amplitude vs. Standardized Frequency',...
            ['Model = ', bubble_model, ', Bubble Radius = ', num2str(R00(rad_ndx)*1e6,'%2.1f'),' \mum'],...
            ['Pulse Magnitude = ',num2str(pacs(pac_ndx)/1e3),' kPa'],...
            ['Reference Rosonance Frequency f_{r} = ',num2str(Finit(rad_ndx,1)/1e6,'%2.2f'),' MHz']});
    end
    xlim([povs(1)/1e3,povs(end)/1e3]);
end

%% 显示不同谐波与环境压力的关系
if length(rads)>=1 && length(povs)>1 && length(frqs)==1 && length(pacs)==1
    fig10130=figure(10130); fig10130.Position = [100 100 1000 800];
    subplot(2,2,1)
    X = povs(:) / 1e3; Y = xsh(:,:,frq_ndx,pac_ndx)';
    plt=plot(X,Y,'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Overpressure (kPa)'); ylabel('Amplitude (dB)');
    if pulse_negative==1
        title({'Subharmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = -',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    else
        title({'Subharmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = ',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    end

    subplot(2,2,2)
    X = povs(:) / 1e3; Y = xfd(:,:,frq_ndx,pac_ndx)';
    plt=plot(X,Y,'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Overpressure (kPa)'); ylabel('Amplitude (dB)');
    if pulse_negative==1
        title({'Fundamental Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = -',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    else
        title({'Fundamental Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = ',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    end
    
    subplot(2,2,3)
    X = povs(:) / 1e3; Y = xuh(:,:,frq_ndx,pac_ndx)';
    plt=plot(X,Y,'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Overpressure (kPa)'); ylabel('Amplitude (dB)');
    if pulse_negative==1
        title({'Ultraharmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = -',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    else
        title({'Ultraharmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = ',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    end
    
    subplot(2,2,4)
    X = povs(:) / 1e3; Y = xhm(:,:,frq_ndx,pac_ndx)';
    plt=plot(X,Y,'Marker','+'); grid on; hold on;
    for k=1:length(plt)
        plt(k).Color=colors{k}; plt(k).LineStyle=styles{k}; plt(k).Marker=mrkers{k}; plt(k).LineWidth=widths{k};
    end
    legend(compose('R0 = %2.1f \x03BCm, \x0394SH = %2.2f dB',rads(:)*1e6, (Y(end,:)-Y(1,:))'),'Location','best');
    xlabel('Overpressure (kPa)'); ylabel('Amplitude (dB)');
    if pulse_negative==1
        title({'Second-Harmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = -',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    else
        title({'Second-Harmonic Amplitude vs. Ambient Overpressure',...
            ['Model = ', bubble_model],...
            ['Pulse Magnitude = ',num2str(Pac/1e3),' kPa, ', 'Pulse Frequency = ',num2str(Frq/1e6),' MHz']});
    end
    
    
end