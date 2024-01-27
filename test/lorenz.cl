void getRHS(const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]) {

    /* State variables */
    realtype x = x_[0];
    realtype y = x_[1];
    realtype z = x_[2];

    /* Parameters */
    realtype r = p_[0];
    realtype s = p_[1];
    realtype b = p_[2];

    /* Differential equations */
    realtype dx = s*(y - x);
    realtype dy = r*x - y - x*z;
    realtype dz = x*y - b*z;

    /* Differential outputs */
    dx_[0] = dx;
    dx_[1] = dy;
    dx_[2] = dz;

    aux_[0] = dx;
}