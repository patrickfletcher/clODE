// A collection of problems from
// https://www.osti.gov/biblio/6111421

void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype m = par[0];
    realtype w = par[1];
    realtype k = par[2];
    realtype H = par[3];

    realtype y1 = var[0];
    realtype y2 = var[1];

    realtype dy1 = y2;
    realtype dy2 = (w - k * y2) / m;

    derivatives[0] = dy1;
    derivatives[1] = dy2;
    aux[0] = y1 - H;
}