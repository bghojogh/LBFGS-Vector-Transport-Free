function [ds]=load_data(dataset_name)
    switch dataset_name
        case 'fisheriris'
            DS_name='fisheriris';
            load('fisheriris');
            XX=meas(1:end,:)';
            DS_labels=species;
            labelsxx(1,1:50)=1;
            labelsxx(51:100)=2;
            labelsxx(101:150)=3;
            XX(end+1,:)=labelsxx;
            XX = XX(:,randperm(size(XX,2)));
            cut=120;
            xTr=XX(1:end-1,1:cut)';
            yTr=XX(end,1:cut)';
            xTe=XX(1:end-1,:)';
            yTe=XX(end,:)';
        case 'banana'
            DS_name='banana';
            DS1_name='banana';
            DS2_name='banana';
            load banana.mat;
            X=banana;
            X = X(:,randperm(size(X,2)));
            labels=X(:,end)';
            X=X(:,1:2)';
            [dx,n] = size(X);
            cut=dx*2/3;
            xTr=X(1:end-1,1:cut)';
            yTr=X(end,1:cut)';
            xTe=X(1:end-1,:)';
            yTe=X(end,:)';
        case 'mnist'
            DS_name='mnist';
            load mnist.mat;
            xTe=XTe';
            xTr=XTr';
            yTr=yTr';
            yTe=yTe';
        case 'usps'
            %http://www.cad.zju.edu.cn/home/dengcai/Data/MLData.html
            %http://www.cad.zju.edu.cn/home/dengcai/Data/Examples.html
            DS_name='usps';
            load('usps');
            xTr=fea(1:7291,:);
            xTe=fea(7292:end,:);
            yTr=gnd(1:7291,:);
            yTe=gnd(7292:end,:);
            data=load('USPS');
            
%             data_te=load('USPS.t');
            %save('d:\uspsx','fea','gnd');
%             xTr=data.XTr';
%             yTr=data.yTr';
%             xTe=data.XTe';
%             yTe=data.yTe';
     
%             DS_name='USPS.mat';
%             load('USPS.mat');
%             xTr=XTr';
%             xTe=XTe';
%             yTr=yTr';
%             yTe=yTe';
        case 'vehicle'
            DS_name='vehicle';
            load('vehicle.mat');
        otherwise
            DS_name='vehicle';
            load('vehicle.mat');       
    end
    s1=size(xTr,1);
    s2=size(xTe,1);
    s3=size(yTr,1);
    s4=size(yTe,1);
    if s1>0 && s2>0 && s3>0 && s4>0
        fprintf("\ndata set "+DS_name+" loaded: xTr=%d xTe=%d yTr=%d yTe=%d",s1,s2,s3,s4);
        ds.xTr=xTr;
        ds.yTr=yTr;
        ds.xTe=xTe;
        ds.yTe=yTe;
    else
        fprintf('\ndata set '+DS_name+' could not be load!!!');
    end
end