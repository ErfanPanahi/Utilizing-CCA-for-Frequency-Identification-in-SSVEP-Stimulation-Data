% be name khoda
% CCA

% Generation of templates
Fs=250;
T=size(data,2);
ts=1/Fs;
t = 0:ts:(T-1)*ts;
freq_len=length(freq);
TEMP = cell(1,freq_len);
for i = 1:freq_len
    temp = [sin(2*pi*freq(i)*t); cos(2*pi*freq(i)*t)];
    for k = 2:7
        if (freq(i)*k > 40)
            break;
        end
        temp = [temp; sin(2*pi*k*freq(i)*t); cos(2*pi*k*freq(i)*t)];
    end
    TEMP{i} = temp;
end


% CCA
trial_num=size(data,3);
label_estimation=zeros(size(label));
for k=1:trial_num
    
    Y=data(:,:,k);
    Ry=Y*Y';
    ro=zeros(1,freq_len);
    for i=1:freq_len
        X=TEMP{i};
        Rx=X*X';
        Rxy=X*Y';
        Ryx=Y*X';
        
        SIGMA1=(Rx^-0.5)*Rxy*(Ry^-1)*Ryx*(Rx^-0.5);
        [V,LANDA]=eig(SIGMA1);
        [landa,pos]=sort(diag(LANDA),'descend');
        V=V(:,pos);
        c=V(:,1);
        ro(i)=landa(1);
        
        SIGMA2=(Ry^-0.5)*Ryx*(Rx^-1)*Rxy*(Ry^-0.5);
        [V,LANDA]=eig(SIGMA2);
        [landa,pos]=sort(diag(LANDA),'descend');
        V=V(:,pos);
        d=V(:,1);
        ro(i)=landa(1);
    end
    
    [~,ind]=max(ro);
    label_estimation(k)=freq(ind);
end

