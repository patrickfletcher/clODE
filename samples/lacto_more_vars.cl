
void getRHS(const realtype t, const realtype x_[], const realtype p_[], realtype dx_[], realtype aux_[], const realtype w_[]) {
realtype minf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-20.0)-x_[0])/RCONST(12.0)));
realtype mtinf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-38.0)-x_[0])/RCONST(6.0)));
realtype htinf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-56.0)-x_[0])/RCONST(-5.0)));
realtype ninf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-5.0)-x_[0])/RCONST(10.0)));
realtype sinf=x_[7]*x_[7]/(x_[7]*x_[7]+RCONST(0.40)*RCONST(0.40));
realtype nkinf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-65.0)-x_[0])/RCONST(-8.0)));
realtype binf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-20.0)-x_[0])/RCONST(2.0)));
realtype v_na_inf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-15.0)-x_[0])/RCONST(5.0)));
realtype hnainf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-60.0)-x_[0])/RCONST(-10.0)));
realtype ainf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-20.0)-x_[0])/RCONST(10.0)));
realtype hinf=RCONST(1.0)/(RCONST(1.0)+exp((RCONST(-60.0)-x_[0])/RCONST(-5.0)));
realtype i_cal=p_[0]*x_[2]*(x_[0]-RCONST(60.0));
realtype i_cat=p_[1]*mtinf*x_[5]*(x_[0]-RCONST(60.0));
realtype i_k=p_[2]*x_[1]*(x_[0]-RCONST(-75.0));
realtype i_sk=p_[3]*sinf*(x_[0]-RCONST(-75.0));
realtype i_kir=p_[4]*nkinf*(x_[0]-RCONST(-75.0));
realtype i_bk=p_[5]*x_[3]*(x_[0]-RCONST(-75.0));
realtype i_nav=p_[6]*v_na_inf*v_na_inf*v_na_inf*x_[6]*(x_[0]-RCONST(75.0));
realtype i_a=p_[7]*ainf*x_[4]*(x_[0]-RCONST(-75.0));
realtype i_leak=p_[8]*(x_[0]-p_[10]);
realtype i_noise=p_[18]*w_[0];
realtype i=i_cal+i_cat+i_k+i_sk+i_kir+i_bk+i_nav+i_a+i_leak+i_noise;
dx_[0]=-i/p_[9];
dx_[1]=(ninf-x_[1])/p_[13];
dx_[2]=(minf-x_[2])/p_[11];
dx_[3]=(binf-x_[3])/p_[14];
dx_[4]=(hinf-x_[4])/p_[15];
dx_[5]=(htinf-x_[5])/p_[12];
dx_[6]=(hnainf-x_[6])/p_[16];
dx_[7]=-RCONST(0.010)*(RCONST(0.00150)*i_cal+p_[17]*x_[7]);
aux_[0]=i_noise;
}

