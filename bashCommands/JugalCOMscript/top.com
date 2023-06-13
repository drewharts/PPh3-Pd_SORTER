%nprocshared=12
%mem=5GB
#p m06/def2svp nosymm  freq=noraman opt=(noeigentest) Pop=(Hirshfeld,NPA)

xyz

1 1
