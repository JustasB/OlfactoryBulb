
clc
clear all;
clf;


UUd= 0.8;
DDd= 0.9;
FFd= 0.1;

UUf= 0.1;
DDf= 0.1;
FFf= 0.9;

%pd= 1.;

%--------------------------------------------------------------------------
JJ = 4.;
II = 5.;

RANGE = 0:0.02:1;
RR = [];
for pd=RANGE
   [Fm Qm r] = FullMixedModel_fixed_point(JJ, II, UUd, DDd, FFd, UUf, DDf, FFf, pd);
   RR = [RR, max(r)];
end

pd=RANGE;
figure(2); clf; hold on;
P1 = plot(pd, RR);
set(gca, 'XLim', [0 1], 'box', 'on', 'FontName', 'Arial', 'FontSize', 20)
xlabel('{P_d}', 'FontName', 'Arial', 'FontSize', 20);

%--------------------------------------------------------------------------

JJ = 10.;
II = 0.;
RR = [];
for pd=RANGE,
   [Fm Qm r] = FullMixedModel_fixed_point(JJ, II, UUd, DDd, FFd, UUf, DDf, FFf, pd);
   RR = [RR, max(r)];
end

pd=RANGE;
P2 = plot(pd, RR);
set(gca, 'XLim', [0 1], 'box', 'on', 'FontName', 'Arial', 'FontSize', 20)
xlabel('{P_d}', 'FontName', 'Arial', 'FontSize', 20);

set(gcf, 'Color', [1 1 1]);
set(P1, 'Color', [0 0 0], 'LineWidth', 2);
set(P2, 'Color', [0 0 0]+0.5, 'LineWidth', 2);
% %ylabel('dh/dt', 'FontName', 'Arial', 'FontSize', 20);
print(gcf, '-depsc2', '-loose', 'mixed.eps');
print(gcf, '-dpng', 'mixed.png');
%--------------------------------------------------------------------------


