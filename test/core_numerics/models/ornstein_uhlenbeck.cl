void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype mu = par[0];
    realtype sigma = par[1];
    realtype x = var[0];

    derivatives[0] = (mu - x) + sigma * wiener[0];
}
