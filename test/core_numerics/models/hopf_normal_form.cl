void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype mu = par[0];
    realtype omega = par[1];
    realtype x = var[0];
    realtype y = var[1];
    realtype r2 = x * x + y * y;

    derivatives[0] = mu * x - omega * y - r2 * x;
    derivatives[1] = omega * x + mu * y - r2 * y;

    aux[0] = r2;
}
