void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype mu = par[0];
    realtype x = var[0];
    realtype y = var[1];

    realtype dx = y;
    realtype dy = mu * (1 - x * x) * y - x;

    derivatives[0] = dx;
    derivatives[1] = dy;
}