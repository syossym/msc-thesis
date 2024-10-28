function Plot_Field(Mat_name)

global Grid;

if length(Mat_name) == 1
    Mat_name = Build_Field_Cavity(Mat_name);
end    
    
imagesc(Grid.axis,Grid.axis,abs(Mat_name))
%surf(Grid.axis,Grid.axis,abs(Mat_name))
shading interp
title('Field in Amplitude Profile')
axis tight
axis square
axis xy
view([0 90])