function METHODS = get_methods(DIM, K, NDATA)
% field names of METHODS define the methods used in the experiments

if nargin < 3
    NDATA = DIM^2*100;
end

if K > 1
    
    it = 50;
    if false
        METHODS.SGDf11.legend = 'SGD (it=50)';
        METHODS.SGDf11.info_fields = struct('iter',[],'gradnorm',[],'cost',[],'time',[],'theta',[],'ll',[]);
        METHODS.SGDf11.ComponentD = mvn2factorytmp(DIM);
        METHODS.SGDf11.options.solver = 'sgd';
        METHODS.SGDf11.options.sgd.stepsize = 1;
        METHODS.SGDf11.options.sgd.batchnum = NDATA/DIM;
        METHODS.SGDf11.options.sgd.epoch = it;
        METHODS.SGDf11.options.sgd.base = 10^(-3/METHODS.SGDf11.options.sgd.batchnum/it);
        METHODS.SGDf11.options.sgd.momentum = 0;
        METHODS.SGDf11.options.sgd.svrg = false;
        METHODS.SGDf11.options.penalize = true;
    end
    
    it = 20;
    if false
        METHODS.SGDf9.legend = 'SGD (it=20)';
        METHODS.SGDf9.info_fields = struct('iter',[],'gradnorm',[],'cost',[],'time',[],'theta',[],'ll',[]);
        METHODS.SGDf9.ComponentD = mvn2factorytmp(DIM);
        METHODS.SGDf9.options.solver = 'sgd';
        METHODS.SGDf9.options.sgd.stepsize = 1;
        METHODS.SGDf9.options.sgd.batchnum = NDATA/DIM;
        METHODS.SGDf9.options.sgd.epoch = it;
        METHODS.SGDf9.options.sgd.base = 10^(-3/METHODS.SGDf9.options.sgd.batchnum/it);
        METHODS.SGDf9.options.sgd.momentum = 0;
        METHODS.SGDf9.options.sgd.svrg = false;
        METHODS.SGDf9.options.penalize = true;
    end
    
    it = 5;
    if false
        METHODS.SGDf12.legend = 'SGD (it=5)'; %m=5d,
        METHODS.SGDf12.info_fields = struct('iter',[],'gradnorm',[],'cost',[],'time',[],'theta',[],'ll',[]);
        METHODS.SGDf12.ComponentD = mvn2factorytmp(DIM);
        METHODS.SGDf12.options.solver = 'sgd';
        METHODS.SGDf12.options.sgd.stepsize = 1;
        METHODS.SGDf12.options.sgd.batchnum = NDATA/DIM;
        METHODS.SGDf12.options.sgd.epoch = it;
        METHODS.SGDf12.options.sgd.base = 10^(-3/METHODS.SGDf12.options.sgd.batchnum/it);
        METHODS.SGDf12.options.sgd.momentum = 0;
        METHODS.SGDf12.options.sgd.svrg = false;
        METHODS.SGDf12.options.penalize = true;
    end
    
    
    if false
        METHODS.EM1.legend = 'EM, Usual MVN';
        METHODS.EM1.info_fields = struct('iter',[],'cost',[],'time',[],'theta',[],'ll',[]);
        METHODS.EM1.ComponentD = mvnfactory(DIM);
        METHODS.EM1.options.solver = 'default';
        METHODS.EM1.options.penalize = true;
    end
    
    % Reparametrized EM and LBFGs are a bit slow
    if false
        METHODS.EM2.legend = 'EM, Reparameterized MVN';
        METHODS.EM2.info_fields = struct('iter',[],'cost',[],'time',[],'theta',[],'ll',[]);
        METHODS.EM2.ComponentD = mvn2factory(DIM);
        METHODS.EM2.options.solver = 'default';
        % added by Reza
        METHODS.EM2.options.penalize=false;
    end
    
    if true
        METHODS.LBFGS1.legend = 'LBFGS, Usual MVN VTFree';
        METHODS.LBFGS1.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'theta',[],'ll',[]);
        METHODS.LBFGS1.ComponentD = mvnfactory_VTFree(DIM);
        METHODS.LBFGS1.options.solver = 'lbfgs';
        % added by Reza
        METHODS.LBFGS1.options.penalize=false;
    end
    
    if true
        METHODS.LBFGS2.legend = 'LBFGS, Usual MVN';
        METHODS.LBFGS2.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'theta',[],'ll',[]);
        METHODS.LBFGS2.ComponentD = mvnfactory(DIM);
        METHODS.LBFGS2.options.solver = 'lbfgs';
        % added by Reza
        METHODS.LBFGS2.options.penalize=false;
        
    end

    if false
        METHODS.LBFGS3.legend = 'LBFGS, Reparameterized MVN';
        METHODS.LBFGS3.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'theta',[],'ll',[]);
        METHODS.LBFGS3.ComponentD = mvn2factory(DIM);
        METHODS.LBFGS3.options.penalize = true;
        METHODS.LBFGS3.options.solver = 'lbfgs';
    end
    
    if false
        METHODS.LBFGS3.legend = 'LBFGS, Factorized MVN';
        METHODS.LBFGS3.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'theta',[],'ll',[]);
        METHODS.LBFGS3.ComponentD = mvn2trilfactory(DIM);
        METHODS.LBFGS3.options.solver = 'lbfgs';
    end
    
    if false
        METHODS.CG3.legend = 'CG, Factorized MVN';
        METHODS.CG3.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'beta',[],'theta',[],'ll',[]);
        METHODS.CG3.ComponentD = mvn2trilfactory(DIM);
        METHODS.CG3.options.solver = 'cg';
    end
    
    
    if false
        METHODS.CG2.legend = 'CG, Reparameterized MVN';
        METHODS.CG2.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'beta',[],'theta',[],'ll',[]);
        METHODS.CG2.ComponentD = mvn2factory(DIM);
        METHODS.CG2.options.penalize = true;
        METHODS.CG2.options.solver = 'cg';
    end
    
    if false
        METHODS.CG1.legend = 'CG, Usual MVN';
        METHODS.CG1.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'beta',[],'theta',[],'ll',[]);
        METHODS.CG1.ComponentD = mvnfactory(DIM);
        METHODS.CG1.options.solver = 'cg';
    end
    
    return
else
    
    METHODS.CG.legend = 'Reparameterized MVN';
    METHODS.CG.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'beta',[],'theta',[],'ll',[]);
    METHODS.CG.ComponentD = mvnfactory2(DIM);
    METHODS.CG.options.solver = 'cg';
    
    METHODS.CG2.legend = 'Usual MVN';
    METHODS.CG2.info_fields = struct('iter',[],'cost',[],'gradnorm',[],'stepsize',[],'time',[],'linesearch',[],'beta',[],'theta',[],'ll',[]);
    METHODS.CG2.ComponentD = mvnfactory(DIM);
    METHODS.CG2.options.solver = 'cg';
    
end

