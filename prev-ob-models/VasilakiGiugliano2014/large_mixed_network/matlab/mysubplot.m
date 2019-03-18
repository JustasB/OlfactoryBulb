function HH = mysubplot(varargin)
%
% HH = mysubplot(N,M,D,Dx,Dy)
%

clf;

if (nargin == 2)
 D = 0.01;
 Dx = 0.01;       Dy = 0.01;
 N = varargin{1};
 M = varargin{2};
else
 N  = varargin{1};
 M  = varargin{2};
 D  = varargin{3};
 Dx = varargin{4};
 Dy = varargin{5};
end %

HH = [];

XMIN = 0 + D;    YMIN = 0 + D;
XMAX = 1 - D;    YMAX = 1 - D;
Lx   = ((XMAX - XMIN) - (M-1)*Dx)/ M;
Ly   = ((YMAX - YMIN) - (N-1)*Dy)/ N;

K = N*M;
for i=1:K,

 x = mod(i-1,M);     X = XMIN + x * (Lx + Dx);
 y = fix((i-1)/M);   Y = YMAX - (y+1) * Ly - y * Dy;

 %H = subplot(N,M,i);
 HH(i) = axes;
 set(HH(i),'Position',[X, Y, Lx, Ly]);
% text(0.5,0.5,num2str(i),'FontSize',8); pause;
end
