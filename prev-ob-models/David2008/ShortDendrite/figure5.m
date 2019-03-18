function figure5(n_trial)

% n_trial :     number of different noise seed.

close all
doublim = 2;
n_fig=15;

for s=1:3
    if s==1
        pa = dlmread('fig5_0.500000.dat');
    elseif s==2
        clear('pa')
        pa = dlmread('fig5_0.800000.dat');
    elseif s==3
        clear('pa')
        pa = dlmread('fig5_1.800000.dat'); 
    end
    
    for i=1:11
        isi{i}=[];
    end    
    for i=1:n_trial
        for j=1:11
            isi{j}=[isi{j} diff(pa(j+11*(i-1),:))];
        end
    end
    isisel{1}=isi{1}(find(isi{1}>doublim));
    isisel{2}=isi{2}(find(isi{2}>doublim+2));
    isisel{3}=isi{3}(find(isi{3}>doublim+3));
    isisel{4}=isi{4}(find(isi{4}>doublim+6));
    isisel{5}=isi{5}(find(isi{5}>doublim+10));
    isisel{6}=isi{6}(find(isi{6}>doublim+12));
    isisel{7}=isi{7}(find(isi{7}>doublim+20));
    isisel{8}=isi{8}(find(isi{8}>doublim+25));
    isisel{9}=isi{9}(find(isi{9}>doublim+30));
    isisel{10}=isi{10}(find(isi{10}>doublim+32));
    isisel{11}=isi{11}(find(isi{11}>doublim+35));
%     
%     isisel{1}=isi{1}(find(isi{1}>doublim+0));
%     isisel{2}=isi{2}(find(isi{2}>doublim+0));
%     isisel{3}=isi{3}(find(isi{3}>doublim+0));
%     isisel{4}=isi{4}(find(isi{4}>doublim+0));
%     isisel{5}=isi{5}(find(isi{5}>doublim+0));
%     isisel{6}=isi{6}(find(isi{6}>doublim+0));
%     isisel{7}=isi{7}(find(isi{7}>doublim+0));
%     isisel{8}=isi{8}(find(isi{8}>doublim+0));
%     isisel{9}=isi{9}(find(isi{9}>doublim+0));
%     isisel{10}=isi{10}(find(isi{10}>doublim+0));
%     isisel{11}=isi{11}(find(isi{11}>doublim+0));
    
    %if s==2
        % automatic figure creation
        for i=[1,4]
            fit_ML_normal(isisel{i},[0:250]);hold on;
            a = histc(isisel{i},[0:250],1);a = a/sum(a);
            bar([0:250],a,1,'stack','FaceColor',[0.7,0.7,0.7]);
            axis([0 100 0 0.06]);set(gca,'FontSize',12);set(gca,'YTick',[0,0.03,0.06]);
            xlabel('ISI (ms)','FontSize',12)
        end
    %end
    for i=1:11
        delai(i) = mean(isisel{i});
        stdd(i) = std(isisel{i});
    end
    figure(n_fig);hold on
    plot(0.02 * (10.^(0.3775*[0:9])),[delai([1:10])],'ko-')

    figure(n_fig+1);hold on
    plot(0.02 * (10.^(0.3775*[0:9])),[stdd([1:10])],'ko-')
end

figure(n_fig);hold on
set(gca,'XScale','log')
set(gca,'XTick',[0 0.05 0.1 0.5 2 5 10 20 50])
set(gca,'XTickLabel',[0 0.05 0.1 0.5 2 5 10 20 50])
xlabel('g_{GABA}')
ylabel('\mu_{ISI} (ms)','Fontsize',14)
axis([0 50 20 80])
    
figure(n_fig+1);hold on
set(gca,'XScale','log')
set(gca,'XTick',[0 0.05 0.1 0.5 2 5 10 20 50])
set(gca,'XTickLabel',[0 0.05 0.1 0.5 2 5 10 20 50])
xlabel('g_{GABA}')
ylabel('\sigma_{ISI} (ms)','Fontsize',14)
axis([0 50 6.5 25])
