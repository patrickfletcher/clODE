void getRHS(const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]) {

    /* State variables */
    realtype x = x_[0];
    realtype y = x_[1];

    /* Parameters */
    realtype mu = p_[0];

    /* Noise terms */

    /* Core equations */
    /*  This comment will be preserved; */

    /* Auxiliary equations */
    realtype x2 = x*x;
    realtype x3 = x*x*x;
    realtype y4 = y*y*y*y;
    realtype x5 = pown(x, 5);
    realtype negx = pow(x, -2);
    realtype negy = pow(y, -1.5f);

    /* Differential equations */
    realtype dy = mu * (1 - x*x) * y - x;
    realtype dx = y;

    /* Auxiliary outputs */
    aux_[0] = x2;
    aux_[1] = x3;
    aux_[2] = y4;
    aux_[3] = x5;
    aux_[4] = negx;
    aux_[5] = negy;

    /* Differential outputs */
    dx_[0] = dx;
    dx_[1] = dy;
}