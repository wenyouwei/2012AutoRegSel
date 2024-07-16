function sigma = ImageStdDev(y)
h = [0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];       
h = h(end:-1:1); %% flip is used only to give exactly same result as previous version of the code (which was using filter2)
        
z = conv2(y,h,'same');
z=conv2(z,h','same');
sigma = median(abs(z(:)))/.6745;