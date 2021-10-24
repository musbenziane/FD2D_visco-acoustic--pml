clear; close all; clc;

load sg.mat

maxamp = max(max(shot_g));  % max amplitude calculation for clip
clip = 99;  % clip value
ampclip = (1-clip/100)*maxamp;
data    = shot_g(1:5000,50:end-50);  
data(data > ampclip)= ampclip; % clipping positive
data(data < -ampclip)= -ampclip; % clipping positive
cmap=colormap(redblue);
imagesc(data);