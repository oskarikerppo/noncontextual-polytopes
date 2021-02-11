import numpy as np
import picos


#function Pwin = nonconIneqUQ2_beta()
AS = 4 # Preparation settings
n = 3 # measurement settings
d=2 #binary outcomes
#-------------------------------------------------------------------------------
# contraint vector declaration.
F0 = [] # contraint for positivity of Gamma. 
F1 = [] # constraints for unitarity level 2
F2 = [] # constraints for summation over results of one measurement.
F3 = [] 
F4 = [] # oblivious!
F5 = [] # measurement equivalences
Gamma = [] #List for SDP variables (moment matrices)
O = 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1 #number of operatos
#-------------------------------------------------------------------------------
# sdpvar declaration.
for a in range(AS):
  sdpvar = picos.HermitianVariable("Moment_matrix_{}".format(a), O)
  Gamma.append(sdpvar)
  F0.append(sdpvar >> 0)


#-------------------------------------------------------------------------------
def idx(b, B, U):
  return 2*(d-1) * b + 2*B + U + 2 - 1# function to return the position of the operator corresponding to Bob's setting (b,B) in a Gamma matrix.

def idx22(b, B, U, b1, B1, U1):
  return 2*n*(d-1)+1 + 2*2*n*(d-1)*(d-1)*b + 2*2*n*(d-1)*B+ 2*n*(d-1)*U + 2*(d-1) * b1 +  2*B1 + U1 + 1 - 1

#-------------------------------------------------------------------------------

# General Matrix constraints F1,F2,F3
for a in range(AS):
  for i in range(O):
    F1.append(Gamma[a][i*O + i] == 1)



#for a = 1:AS
for a in range(AS):
#  for b = 0:(n-1)
  for b in range(n):
#    for B = 0:(d-2)
    for B in range(d-1):
#        for U = 0:1
      for U in range(0, 2):
#            for b1 = 0:(n-1)
        for b1 in range(0, n):
#                for B1 = 0:(d-2)
          for B1 in range(0, d-1):
#                    for U1 = 0:1
            for U1 in range(0, 2):                      
#            for br = 0:(n-1)
              for br in range(n):
#                for Br = 0:(d-2)
                for Br in range(d-1):
#                    for Ur = 0:1
                  for Ur in range(2):
                    if (b == b1) & (B == B1) & (U!=U1):
                      #F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br,Br,Ur),1)];
                      F2.append(Gamma[a][O*idx(br,Br,Ur) + idx22(b,B,U,b1,B1,U1)] == Gamma[a][O*idx(br,Br,Ur)])
                      #F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(1,idx(br,Br,Ur))];
                      F2.append(Gamma[a][O*idx22(b,B,U,b1,B1,U1) + idx(br,Br,Ur)] == Gamma[a][idx(br,Br,Ur)])
                    elif (b == br) & (B==Br) & (U==Ur):

                      #F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(1,idx(b1,B1,U1))];
                      #F2.append(Gamma[a][O*idx(br,Br,Ur) + idx22(b,B,U,b1,B1,U1)] == Gamma[a][idx(b1,B1,U1)])
                      #F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(idx(b1,B1,U1),1)];
                      #F2.append(Gamma[a][O*idx22(b,B,U,b1,B1,U1) + idx(br,Br,Ur)] == Gamma[a][O*idx(b1,B1,U1)])
                      pass
                    elif (b1 == br) & (B1==Br) & (U1==Ur):
                      #F2=[F2;Gamma{a}(idx(br,Br,Ur),idx22(b,B,U,b1,B1,U1))==Gamma{a}(1,idx(b,B,U))];
                      F2.append(Gamma[a][O*idx(br,Br,Ur) + idx22(b,B,U,b1,B1,U1)] == Gamma[a][idx(b,B,U)])
                      #F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx(br,Br,Ur))==Gamma{a}(idx(b,B,U),1)];
                      F2.append(Gamma[a][O*idx22(b,B,U,b1,B1,U1) + idx(br,Br,Ur)] == Gamma[a][O*idx(b,B,U)])
                    #end
                  #for br1 = 0:(n-1)
                  for br1 in range(n):
                    #for Br1 = 0:(d-2)
                    for Br1 in range(0,d-1):
                      #for Ur1 = 0:1
                      for Ur1 in range(2):
                        if (br==b) & (Br == B) & (U==Ur):
                            #F2=[F2;Gamma{a}(idx22(br,Br,Ur,br1,Br1,Ur1),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br1,Br1,Ur1),idx(b1,B1,U1))];
                            F2.append(Gamma[a][O*idx22(br,Br,Ur,br1,Br1,Ur1) + idx22(b,B,U,b1,B1,U1)] == Gamma[a][O*idx(br1,Br1,Ur1) + idx(b1,B1,U1)])
                            #F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx22(br,Br,Ur,br1,Br1,Ur1))==Gamma{a}(idx(b1,B1,U1),idx(br1,Br1,Ur1))];
                            F2.append(Gamma[a][O*idx22(b,B,U,b1,B1,U1) + idx22(br,Br,Ur,br1,Br1,Ur1)] == Gamma[a][O*idx(b1,B1,U1) + idx(br1,Br1,Ur1)])
                        elif (br1==b1) & (Br1 == B1) & (Ur1==U1):
                            #F2=[F2;Gamma{a}(idx22(br,Br,Ur,br1,Br1,Ur1),idx22(b,B,U,b1,B1,U1))==Gamma{a}(idx(br,Br,Ur),idx(b,B,U))];
                            F2.append(Gamma[a][O*idx22(br,Br,Ur,br1,Br1,Ur1) + idx22(b,B,U,b1,B1,U1)] == Gamma[a][O*idx(br,Br,Ur) + idx(b,B,U)])
                            #F2=[F2;Gamma{a}(idx22(b,B,U,b1,B1,U1),idx22(br,Br,Ur,br1,Br1,Ur1))==Gamma{a}(idx(b,B,U),idx(br,Br,Ur))];
                            F2.append(Gamma[a][O*idx22(b,B,U,b1,B1,U1) + idx22(br,Br,Ur,br1,Br1,Ur1)] == Gamma[a][O*idx(b,B,U) + idx(br,Br,Ur)])


#-------------------------------------------------------------------------------     
#F4 = [F4;Gamma{1}+Gamma{2}==Gamma{3}+Gamma{4};Gamma{1}+Gamma{2}==Gamma{5}+Gamma{6}];
F4.append(Gamma[0] + Gamma[1] == Gamma[2] + Gamma[3])
#F4.append(Gamma[0] + Gamma[1] == Gamma[4] + Gamma[5])
#-------------------------------------------------------------------------------     

#for a = 1:AS
for a in range(AS):
  #for b = 1: 4*n*n*(d-1)*(d-1) + 2*n*(d-1)+1
  for b in range(O):
    ssum = Gamma[a][O*b + idx(0, 0, 0)] + Gamma[a][O*b + idx(0, 0, 1)] + Gamma[a][O*b + idx(1, 0, 0)] + Gamma[a][O*b + idx(1, 0, 1)] + Gamma[a][O*b + idx(2, 0, 0)] + Gamma[a][O*b + idx(2, 0, 1)];
    tum = Gamma[a][O*idx(0, 0, 0) + b] + Gamma[a][O*idx(0, 0, 1) + b] + Gamma[a][O*idx(1, 0, 0) + b] + Gamma[a][O*idx(1, 0, 1) + b] + Gamma[a][O*idx(2, 0, 0) + b] + Gamma[a][O*idx(2, 0, 1) + b];
    #F5 = [F5;sum==0;tum==0];
    F5.append(ssum == 0)
    F5.append(tum == 0)
    

#----
#-------------------------------------------------------------------------------     
# Making a probability cell.
#Pcell = cell(AS,n,d); # a0, a1, b, B
#for a = 1:AS
#  for b = 0:(n-1)
#    for B = 0:(d-2)
#            Pcell{a,b+1,B+1 } =0.5 + (Gamma{a}(idx(b,B,0),1)+Gamma{a}(idx(b,B,1),1))/4;
         #   F4 = [F4 ; 0.5 + (Gamma{a}(1, idx(b,B,0))+Gamma{a}(1, idx(b,B,1)))/4>=0];
            
#    end
#  end
#end
a_l = []
for a in range(AS):
  b_l = []
  for b in range(n):
    B_l= []
    for B in range(d-1):
      B_l.append(0.5 + (Gamma[a][O*idx(b,B,0)] + Gamma[a][O*idx(b,B,1)])/4)
    b_l.append(B_l)
  a_l.append(b_l)
Pcell = a_l

#-------------------------------------------------------------------------------  
Pwin = [];
#-------------------------------------------------------------------------------  
#Pwin1 = real(Pcell{1,1,1} + Pcell{3,2,1} + Pcell{5,3,1}) ;
Pwin1 = (Pcell[0][0][0] + Pcell[2][1][0] + Pcell[3][2][0]).real
P = picos.Problem()
P.set_objective("max", Pwin1)
P.add_list_of_constraints(F0)
P.add_list_of_constraints(F1)
P.add_list_of_constraints(F2)
P.add_list_of_constraints(F3)
P.add_list_of_constraints(F4)
P.add_list_of_constraints(F5)
P.solve(solver = "cvxopt")
print(Pwin1.value)
Pwin.append(Pwin1)

'''
#-------------------------------------------------------------------------------  
Pwin2 = real(Pcell{1,1,1} + Pcell{2,2,1} + Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin2, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin2)];
#------------------------------------------------------------------------------- 
Pwin3 = real(Pcell{1,1,1} - Pcell{3,1,1} -2 * Pcell{5,1,1} -2 * Pcell{2,2,1} + 2 * Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin3, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin3)];
#------------------------------------------------------------------------------- 
Pwin4 = real(2* Pcell{1,1,1} - Pcell{2,2,1} +2* Pcell{3,2,1}); 
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin4, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin4)];
#------------------------------------------------------------------------------- 
Pwin5 = real(Pcell{1,1,1} - Pcell{5,1,1} +  Pcell{2,2,1} + Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin5, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin5)];
#------------------------------------------------------------------------------- 
Pwin6 = real(Pcell{1,1,1} - Pcell{5,1,1} +  2*Pcell{2,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin6, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin6)];
#-------------------------------------------------------------------------------  
Pwin7 =real(Pcell{1,1,1} - Pcell{4,1,1} -2 * Pcell{5,1,1} -2 * Pcell{2,2,1} + 2 * Pcell{3,2,1} + 2 * Pcell{5,3,1});
diagnostics = optimize([F0;F1;F2;F3;F4;F5;], -Pwin7, sdpsettings('solver', 'sdpt3'))
Pwin = [Pwin;value(Pwin7)];
'''