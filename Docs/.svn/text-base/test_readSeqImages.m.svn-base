% By KS 6/10.


h19=figure(19); 
cla;
colormap(gray)
axis image

%% animate image sequence; full frame
h19=figure(19); 
cla;
colormap(gray)

frames=readSeqImages('090508-23b_0063.seq', 'lastFrame', 100);


hi=imagesc(squeeze(frames(1, :, :)));
axis image
for i=2:size(frames, 1)
   set(hi, 'CData',  squeeze(frames(i, :, :)))
   xlabel(strcat('frame-', num2str(i), '/', num2str(size(frames, 1))));
   fprintf( 1, '\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%6d', i, size(frames, 1));
   %set( gca, 'visible', 'off' )
	pause( .1) % allow plotting to finish
end

%% animate image sequence; roi
h19=figure(19); 
cla;
colormap(gray)

frames=readSeqImages('090508-23b_0063.seq', 'roi', [51 51 100 100], 'lastFrame', 100);


hi=imagesc(squeeze(frames(1, :, :)));
axis image
for i=2:size(frames, 1)
   set(hi, 'CData',  squeeze(frames(i, :, :)))
   xlabel(strcat('frame-', num2str(i), '/', num2str(size(frames, 1))));
   fprintf( 1, '\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%6d', i, size(frames, 1));
   %set( gca, 'visible', 'off' )
	pause( .1) % allow plotting to finish
end



%% linescan
h19=figure(19); 
cla;
colormap(gray)

frames=readSeqImages('090508-23b_0063.seq', 'roi', [1 100 200 1], ...
    'lastFrame', 350);

imagesc(frames)

