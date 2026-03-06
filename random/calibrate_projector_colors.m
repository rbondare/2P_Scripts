parent_dir='K:\fs3-joeschgrp\projector spectra\';%path to grpdrive\projector spectra
correct_lens=true; %correct for lens transmission
norm_by_area=false; %normalize templates by area 
repeated_measure=true; % repeat for newly generated spectrum files

% figure

S_max=375; %max peak S
M_max=508; %max peak M

base_spectrum_file=[parent_dir 'sky_cloud_190113.mat'];
%base_spectrum_file=[parent_dir 'cloud.mat'];
lens_file=[parent_dir 'lens_transmission.mat'];
projector_spectrum_file=[parent_dir 'test.mat'];

load(lens_file,'lens_template');
load(base_spectrum_file)
if correct_lens
    base_spectrum=spectrum.*lens_template;
else
    base_spectrum=spectrum;
end
base_spectrum=base_spectrum./sum(base_spectrum);

nm=wavelength_nmAir;

S_template=calculate_template(S_max,nm);
M_template=calculate_template(M_max,nm);
if norm_by_area
    S_template=S_template./sum(S_template);
    M_template=M_template./sum(M_template);
end
%
base_relativeMS=sum(base_spectrum.*M_template)/sum(base_spectrum.*S_template);
X=categorical({'base','projector'});
b=1;


while b~=3
    projector_spectrum=load(projector_spectrum_file,'spectrum');
    if correct_lens
    projector_spectrum=projector_spectrum.spectrum.*lens_template;
    end
    projector_spectrum=projector_spectrum/sum(projector_spectrum);
        
    projector_relativeMS=sum(projector_spectrum.*M_template)/sum(projector_spectrum.*S_template);
    
    bar(X,[base_relativeMS,projector_relativeMS]);
    title({sprintf('M/S base %.2f, proj %.2f',base_relativeMS,projector_relativeMS);'click left to load new, right for done'})
    ylabel('relative M over S activation')
    if repeated_measure;[~,~,b]=ginput(1);else b=3;end
end
%close
%% plot norm spectra
figure
plot(nm,base_spectrum,nm,projector_spectrum,nm,M_template/sum(M_template),nm,S_template/sum(S_template),nm,lens_template./sum(lens_template));
legend('base','proj','norm M','normS','lens')

% plot channel spectra
figure
plot(nm,base_spectrum.*M_template,nm,base_spectrum.*S_template,nm,projector_spectrum.*M_template,nm,projector_spectrum.*S_template);
legend('baseM','baseS','projM','projS')

%%

function [template] = calculate_template(nm,range)
% nm: enter the peak of the absorption spectrum
% range: vector of an arbitrary range in nm


lambda_max = nm;
x = lambda_max./(range);
x2 =(range);
%% alpha band

% set parameters from Govardoskii et al. 2000
A= 69.7;
alpha=0.8795 + 0.0459*exp(-((lambda_max-300)^2/11940));
B=28;
beta=0.922;
C=-14.9;
gamma=1.104;
D=0.674;

S = 1./(exp(A*(alpha-x))+exp(B*(beta-x))+exp(C*(gamma-x))+D);

%% beta Band

% set parameters from Govardoskii et al. 2000
lambda_b = 189+0.315*lambda_max;
A_b= 0.26;
beta_b=-40.5+0.195*lambda_max;


Sb= A_b*exp(-(((x2-lambda_b)/beta_b).^2));


template = (Sb+S);%./max(Sb+S);
end
