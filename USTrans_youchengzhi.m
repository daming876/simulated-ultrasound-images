clc
clear
close all

filepath=pwd;           %保存当前工作目录
n_pw = (1:25);
Na = n_pw;
alpha_max = deg2rad(16);  % alpha_max = atan(1/2/F_number);
alpha=linspace(-alpha_max,alpha_max,length(n_pw));    % vector of angles [rad]

for i=33:68
picname=['ultrasval',num2str(str2num(int2str(i)),'%05d'),'.bmp'];
I = imread(picname);
I = rgb2gray(I);
[nl,nc] = size(I);
% imshow(I,'InitialMagnification','fit')
% Choose a cardiac phased array
param = getparam('L11-5v');
param.fc = randi([50,130],1)*1e5;
param.fs = 4*param.fc;
param.kerf = randi([2,5],1)*1e-5;
param.width = randi([25,30],1)*1e-5;
param.pitch = param.kerf + param.width;
array_lenth = param.Nelements*param.pitch;
% Create a 2-D distribution of scatterers whose depth is 15 cm. RC will contain the reflection coefficients.
[xs,~,zs,RC] = genscat([NaN 5e-2],param,I);
scale = (max(xs)-min(xs))/array_lenth;
xs=xs/scale;
zs=20/1000 + zs/scale;
% Create a 256x256 polar grid with IMPOLGRID.
xi = linspace(min(xs),max(xs),nc); % in m
zi = linspace(min(zs),max(zs),nl); % in m
[xi,zi] = meshgrid(xi,zi);

IQc = zeros(size(xi),'like',1i); % will contain the compound I/Q
opt.ElementSplitting = 1; % to make simulations faster
opt.WaitBar = false; % no progress bar for SIMUS
param.fs = 4*param.fc; % sampling frequency
param.fnumber = [];
% % figure(2)
% % Display the scatterers
% scatter(xs*1e2,zs*1e2,2,RC,'filled');
% colormap([hot; 1-hot])
% axis off
% set(gca,'XColor','none','box','off')
% % set(gca,'position',[0 0 1 1]);
% % title([int2str(numel(RC)) ' scatterers'])
% % hold on
h = waitbar(0,'SIMUS & DAS...');
for k = 1:length(n_pw)
 dels = txdelay(param,alpha(k)); % transmit delays
 RF = simus(xs,zs,RC,dels,param,opt); % RF simulation
 IQ = rf2iq(RF,param); % I/Q demodulation
 IQb = das(IQ,xi,zi,dels,param); % DAS beamforming
 IQc = IQc+IQb; % compounding
 waitbar(k/length(alpha),h,...
 ['SIMUS & DAS: ' int2str(k) ' of ' int2str(length(n_pw)) ' completed'])
end
close(h)

B = bmode(IQc,60); % log-compressed image
figure
imagesc(xi(1,:)*1e2,zi(:,1)*1e2,B)
set(gcf,'Position',[200 200 nc nl]);%消除白边
set(gca,'Position',[0 0 1 1]);%消除白边
% shading interp, axis equal ij tight
colormap gray
% axis equal tight
box off
axis off

cd('D:\MUST (1)\US VALIDATE(1)')          %把当前工作目录切换到指定文件夹
usimagename=[num2str(str2num(int2str(i)),'%05d'),'.jpg'];
IR = frame2im(getframe(gcf)); %Convert plot to image (true color RGB matrix).
IR = rgb2gray(IR);
J = imresize(IR,[nl,nc], 'bicubic'); %Resize image to resolution 320x240
imwrite(J,usimagename);
save(num2str(str2num(int2str(i)),'%05d'),'param')
cd(filepath)
end