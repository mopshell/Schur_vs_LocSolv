classdef oppDFLL < oppSpot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        mt, model, U, LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = oppDFLL(mt, model, U, LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr)
           m            = model.nrec*model.nsurf*length(model.freq);
           n            = numel(mt);
           op           = op@oppSpot('oppDFLL', m, n);
           op.cflag     = 1;  
           op.linear    = 1;
           op.children  = []; 
           op.sweepflag = 0;
           op.mt        = mt;
           op.U         = U;
           op.model     = model;
           op.LD        = LD;
           op.Gplus1    = Gplus1;
           op.Gdc       = Gdc;
           op.U0dc      = U0dc;
           op.Gplusr    = Gplusr;
           op.Gdcr      = Gdcr;
       end 
       
    end
    
    
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           if mode == 1
                y = DFLL(op.mt, x, op.model, op.U, op.LD, op.Gplus1, op.Gdc, op.U0dc, op.Gplusr, op.Gdcr, 1);
                
           else %adjoint
                y = DFLL(op.mt, x, op.model, op.U, op.LD, op.Gplus1, op.Gdc, op.U0dc, op.Gplusr, op.Gdcr, -1);  
           end
       end %multiply
       
    end %protected methods
    
end %classdef

    
    
    
