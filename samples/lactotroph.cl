
void getRHS(realtype t, realtype x_[], realtype p_[], realtype dx_[], realtype aux_[], realtype w_[]) {
realtype c2=x_[3]*x_[3];
realtype vkdrive=x_[0]-p_[7];
realtype minf=(1.0)/((1.0)+exp(((-20.0)-x_[0])/(12.0)));
realtype ninf=(1.0)/((1.0)+exp(((-5.0)-x_[0])/(10.0)));
realtype ica=p_[0]*minf*(x_[0]-p_[6]);
realtype isk=p_[1]*c2/(c2+(0.40)*(0.40))*vkdrive;
realtype ibk=p_[2]*x_[2]*vkdrive;
realtype ik=p_[3]*x_[1]*vkdrive;
realtype il=p_[4]*(x_[0]-p_[8]);
realtype itot=ica+isk+ibk+ik+il;
dx_[0]=-itot/p_[9];
dx_[1]=(ninf-x_[1])/p_[10];
dx_[2]=((1.0)/((1.0)+exp(((-20.0)-x_[0])/(2.0)))-x_[2])/p_[11];
dx_[3]=-(0.010)*((0.00150)*ica+p_[5]*x_[3]);
aux_[0]=ica;
}
