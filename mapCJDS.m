load('CJDS_dem_LR.mat')
load('CJDS_linked.mat')
figure(1);
ax = gca;

ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;


pcolor(LonLR,LatLR,DEMlr); shading flat; axis equal; colormap gray
hold on;
set(0,'defaultAxesFontName', 'Helvetica')
set(0,'defaultTextFontName', 'Helvetica')

myColorMap=hsv(90);
for i=1:length(Lat1)
  theIntensity = ceil((AvgTheta(i)+90)/2); % Note row,column = y,x, not x,y.
  % figure out the color from that intensity.
  theLineColor = myColorMap(theIntensity, :);
  line( [Lon1(i), Lon2(i)],[Lat1(i), Lat2(i)], 'Color', theLineColor)

end

%colormap(myColorMap)

xlim([-118, -116.5])
ylabel('Latitude')
xlabel('Longitude')
a = colorbar('eastoutside');
%caxis([-90, 90])
%a.Label.String = 'Theta (deg)';
ylim([44.00, 46.5])
