function test_suite = test_oppBlockDiag
%test_oppBlockDiag  Unit tests for the opBlockDiag operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_oppBlockDiag_builtin
%%
   n = randi(10); m = randi(10);
   A = opMatrix(randn(m,m));
   B = opMatrix(randn(n,n));
   D = oppBlockDiag(A,B);
   utest(D,1);
end

function test_oppBlockDiag_prod
%%
   n = randi(10); m = randi(10);
   A = opMatrix(randn(m,m));
   B = opMatrix(randn(n,n));
   D = oppBlockDiag(A,B,1);
   x = drandn(D,2);
   x2 = gather(x);
   assertElementsAlmostEqual( norm([A*x2(1:m,:); B*x2(m+1:end,:)]- D*x),0 )
   assertElementsAlmostEqual( norm([A'*x2(1:m,:); B'*x2(m+1:end,:)]- D'*x),0 )
end

function test_oppBlockDiag_repeat
%%
    n  = randi([2,10]);
    A  = opMatrix(randn(n,n));
    D  = oppBlockDiag(3,A,1);
    E  = opBlockDiag(3,A);
    x  = drandn(D,2);
    x2 = gather(x);    
    assertElementsAlmostEqual( D*x, E*x2 )
    assertElementsAlmostEqual( D'*x, E'*x2 )
end

function test_oppBlockDiag_weights
%%
   n = randi(10); m = randi(10);
   A = randn(m,m);
   B = randn(n,n);
   D = oppBlockDiag([m n],A,B,1);
   A2 = opMatrix(m*A); B2 = opMatrix(n*B);
   E = opBlockDiag(A2,B2);
   x = drandn(D,2);
   x2 = gather(x);
   assertElementsAlmostEqual(D*x, E*x2);
   assertElementsAlmostEqual(D'*x,E'*x2);
end

function test_oppBlockDiag_divide
%%
    n = randi(10); m = randi(10);
    A1 = opMatrix(randn(m,m));
    A2 = opMatrix(randn(n,n));
    D = oppBlockDiag(A1,A2);
    x = drandn(D,2);
    y = D*x;
    xt = D\y;
    assertElementsAlmostEqual(gather(x),gather(xt));
end