
void getRHS(realtype t, realtype x_[], realtype p_[], realtype dx_[], realtype aux_[], realtype w_[]) {
realtype inoise=p_[2]*w_[0];
realtype c2=x_[3]*x_[3];
realtype minf=(1.0f)/((1.0f)+exp(((-20.0f)-x_[0])/(12.0f)));
realtype ninf=(1.0f)/((1.0f)+exp(((-5.0f)-x_[0])/(10.0f)));
realtype finf=(1.0f)/((1.0f)+exp(((-20.0f)-x_[0])/(2.0f)));
realtype ica=(1.50f)*minf*(x_[0]-(60.0f));
realtype vkdrive=x_[0]-(-75.0f);
realtype isk=(3.0f)*c2/(c2+(0.40f)*(0.40f))*vkdrive;
realtype ibk=p_[0]*x_[2]*vkdrive;
realtype ik=(2.0f)*x_[1]*vkdrive;
realtype il=(0.050f)*(x_[0]-(-50.0f));
realtype itot=ica+isk+ibk+ik+il+inoise;
dx_[0]=-itot/(10.0f);
dx_[1]=(ninf-x_[1])/(30.0f);
dx_[2]=(finf-x_[2])/p_[1];
dx_[3]=-(0.010f)*((0.00150f)*ica+(0.20f)*x_[3]);
aux_[0]=ica;
}

