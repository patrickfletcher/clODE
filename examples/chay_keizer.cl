void getRHS(const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]) {

    /* State variables */
    realtype v = x_[0];
    realtype n = x_[1];
    realtype c = x_[2];

    /* Parameters */
    realtype gca = p_[0];
    realtype gkca = p_[1];
    realtype kpmca = p_[2];
    realtype gk = RCONST(3000.0);

    realtype vca = RCONST(25.0);
    realtype vk = RCONST(-75.0);
    realtype cm = RCONST(5300.0);
    realtype alpha = RCONST(4.5e-6);
    realtype fcyt = RCONST(0.01);
    realtype kd = RCONST(0.4);
    realtype vm = RCONST(-20.0);
    realtype sm = RCONST(12.0);
    realtype vn = RCONST(-16.0);
    realtype sn = RCONST(5.0);
    realtype taun = RCONST(20.0);
    
    /* activation functions */
    realtype minf=RCONST(1.0)/(RCONST(1.0)+exp((vm-v)/sm));
    realtype ninf=RCONST(1.0)/(RCONST(1.0)+exp((vn-v)/sn));
    realtype omega=pown(c, 2)/(pown(c, 2)+pown(kd, 2));

    realtype ica=gca*minf*(v-vca);
    realtype ik=gk*n*(v-vk);
    realtype ikca=gkca*omega*(v-vk);

    /* Differential equations */
    dx_[0] = -(ica + ik + ikca)/cm;
    dx_[1] = (ninf-n)/taun;
    dx_[2] = fcyt*(-alpha*ica - kpmca*c);

}