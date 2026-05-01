void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype a = par[0];
    realtype b = par[1];
    realtype x = var[0];
    realtype y = var[1];

    derivatives[0] = -a * x;
    derivatives[1] = -b * y;
}
