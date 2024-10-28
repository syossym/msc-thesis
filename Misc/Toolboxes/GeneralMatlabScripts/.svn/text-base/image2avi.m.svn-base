function image2avi(List,OutputName,FPS)
%----------------------------------------------------------------
% image2avi function    Create avi file from list of images
% Input  : - File name containing list of images (e.g., in jpg)
%          - Output avi file name.
%          - frame per sec, default is 25.
% Example : image2avi('image.list','out.avi',5);
% Tested : Matlab 6.1
%     By : Eran O. Ofek        May 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%----------------------------------------------------------------

if (nargin==2),
   FPS = 25;
elseif (nargin==3),
   % do nothing
else
   error('Illigal number of input arguments');
end




FID = fopen(List,'r');

%fig=figure;
%set(fig,'DoubleBuffer','on');
MovH = avifile(OutputName,'compression','none','fps',FPS);

%N = wc(List,'l');
for I=1:1:110,
%while (feof(FID)~=1),

    ImageName = fgetl(FID);
    Image     = imread(ImageName);
    H = imshow(uint8(Image));
    set(H,'EraseMode','xor');
    F = getframe(gca);
    MovH = addframe(MovH,F);
end
MovH = close(MovH);

fclose(FID);


