function Pwin = nonconIneqUQ2_beta()
AS = 4;
n = 6;
d=2;
%-------------------------------------------------------------------------------
% contraint vector declaration.
F0 = []; % contraint for positivity of Gamma. 
F1 = []; % constraints for unitarity level 2
F2 = []; % constraints for summation over results of one measurement.
F3 = []; 
F4 = []; % oblivious!
F5 = []; % measurement equivalences
%-------------------------------------------------------------------------------
% sdpvar declaration.
for a = 1:AS
  Gamma{a}=sdpvar( 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1, 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1,'hermitian','real');
  F0=[F0;Gamma{a}>=0];
end
%-------------------------------------------------------------------------------
idx = @(b, B, U) 2*(d-1) * b + 2*B + U + 2 ; % function to return the position of the operator corresponding to Bob's setting (b,B) in a Gamma matrix.
idx22 = @(b, B, U, b1, B1, U1) 2*n*(d-1)+1 + 2*2*n*(d-1)*(d-1)*b + 2*2*n*(d-1)*B+ 2*n*(d-1)*U + 2*(d-1) * b1 +  2*B1 + U1 + 1 
%-------------------------------------------------------------------------------
% General Matrix constraints F1,F2,F3
for a = 1:AS
    for i = 1:4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1
        F1 = [F1;Gamma{a}(i,i)==1];
  end
end
for a = 1:AS
  for b = 0:(n-1)
    for B = 0:(d-2)
        for U = 0:1
            for b1 = 0:(n-1)
                for B1 = 0:(d-2)
                    for U1 = 0:1
                        
            for br = 0:(n-1)
                for Br = 0:(d-2)
                    for Ur = 0:1
                        if (b == b1) && (B == B1) && (U~=U1)
                            F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br,Br,Ur),1)];
                            F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(1,idx(br,Br,Ur))];
                        elseif (b == br) && (B==Br) && (U==Ur)
                       %    F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(1,idx(b1,B1,U1))];
                         %  F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(idx(b1,B1,U1),1)];
                        elseif (b1 == br) && (B1==Br) && (U1==Ur)
                           F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(1,idx(b,B,U))];
                           F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(idx(b,B,U),1)];
                        end
                        for br1 = 0:(n-1)
                for Br1 = 0:(d-2)
                    for Ur1 = 0:1
                        if (br==b) && (Br == B) && (U==Ur)
                            F2=[F2;Gamma{a}(idx22(br,Br,Ur,br1,Br1,Ur1),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br1,Br1,Ur1),idx(b1,B1,U1))];
                            F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx22(br,Br,Ur,br1,Br1,Ur1))==Gamma{a}(idx(b1,B1,U1),idx(br1,Br1,Ur1))];
                        elseif (br1==b1) && (Br1 == B1) && (Ur1==U1)
                            F2=[F2;Gamma{a}(idx22(br,Br,Ur,br1,Br1,Ur1),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br,Br,Ur),idx(b,B,U))];
                            F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx22(br,Br,Ur,br1,Br1,Ur1))==Gamma{a}(idx(b,B,U),idx(br,Br,Ur))];
                        end
                    end
                end
                        end
                end
                end
            end
            
            
                    end
                end
            end
        end
    end
  end
end


%-------------------------------------------------------------------------------     
%F4 = [F4;Gamma{1}+Gamma{2}==Gamma{3}+Gamma{4}];
%-------------------------------------------------------------------------------     
for a = 1:AS
    for b = 1: 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1
               sum = Gamma{a}(b, idx(0, 0, 0)) + Gamma{a}(b, idx(0, 0, 1)) + Gamma{a}(b, idx(1, 0, 0)) +Gamma{a}(b, idx(1, 0, 1)) + Gamma{a}(b, idx(2, 0, 0)) +Gamma{a}(b, idx(2, 0, 1));
               tum = Gamma{a}(idx(0, 0, 0), b) + Gamma{a}(idx(0, 0, 1), b) + Gamma{a}(idx(1, 0, 0), b) +Gamma{a}(idx(1, 0, 1), b) + Gamma{a}(idx(2, 0, 0), b) +Gamma{a}(idx(2, 0, 1), b);
               F5 = [F5;sum==0;tum==0];
    end
end

%----
%-------------------------------------------------------------------------------     
% Making a probability cell.
Pcell = cell(AS,n,d); % a0, a1, b, B
for a = 1:AS
  for b = 0:(n-1)
    for B = 0:(d-2)
            Pcell{a,b+1,B+1 } =0.5 + (Gamma{a}(idx(b,B,0),1)+Gamma{a}(idx(b,B,1),1))/4;
         %   F4 = [F4 ; 0.5 + (Gamma{a}(1, idx(b,B,0))+Gamma{a}(1, idx(b,B,1)))/4>=0];
            
    end
  end
end
%-------------------------------------------------------------------------------  
Pwin = [];
%-------------------------------------------------------------------------------  
Pwin1 = real(-Pcell{1,1,1}-Pcell{1,2,1}-Pcell{1,3,1}+Pcell{2,1,1}-Pcell{2,4,1}-Pcell{2,5,1}+Pcell{3,2,1}+Pcell{3,4,1}-Pcell{3,6,1}+Pcell{4,3,1}+Pcell{4,5,1}+Pcell{4,6,1});

diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin1, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin1)];