
realtype x_inf(const realtype V,
               const realtype V_x,
               const realtype s_x) {
    return RCONST(1.0) / (RCONST(1.0) + exp((V_x - V) / s_x));
}


realtype s_inf(const realtype c,
               const realtype k_s) {
    realtype c2 = c * c;
    return c2 / (c2 + k_s * k_s);
}


void getRHS(const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]) {
    realtype V = x_[0];
    realtype n = x_[1];
    realtype m = x_[2];
    realtype b = x_[3];
    realtype h = x_[4];
    realtype h_T = x_[5];
    realtype h_Na = x_[6];
    realtype c = x_[7];

    realtype g_CaL = p_[0];
    realtype g_CaT = p_[1];
    realtype g_K = p_[2];
    realtype g_SK = p_[3];
    realtype g_Kir = p_[4];
    realtype g_BK = p_[5];
    realtype g_NaV = p_[6];
    realtype g_A = p_[7];
    realtype g_leak = p_[8];
    realtype C_m = p_[9];
    realtype E_leak = p_[10];
    realtype tau_m = p_[11];
    realtype tau_ht = p_[12];
    realtype tau_n = p_[13];
    realtype tau_BK = p_[14];
    realtype tau_h = p_[15];
    realtype tau_hNa = p_[16];
    realtype k_c = p_[17];


    realtype E_Ca = RCONST(60.);
    realtype E_K = RCONST(-75.);
    realtype E_Na = RCONST(75.);

    realtype V_m = RCONST(-20.);
    realtype V_mt = RCONST(-38.);
    realtype V_ht = RCONST(-56.);
    realtype V_n = RCONST(-5.);
    realtype V_k = RCONST(-65.);
    realtype V_b = RCONST(-20.);
    realtype V_a = RCONST(-20.);
    realtype V_h = RCONST(-60.);
    realtype V_mNa = RCONST(-15.);
    realtype V_hNa = RCONST(-60.);

    realtype s_m = RCONST(12.);
    realtype s_mt = RCONST(6.);
    realtype s_ht = RCONST(-5.);
    realtype s_n = RCONST(10.);
    realtype s_k = RCONST(-8.);
    realtype s_b = RCONST(2.);
    realtype s_a = RCONST(10.);
    realtype s_h = RCONST(-5.);
    realtype s_mNa = RCONST(5.);
    realtype s_hNa = RCONST(-10.);

    realtype f_c = RCONST(0.01);
    realtype alpha = RCONST(0.0015);
    realtype k_s = RCONST(0.4);

    realtype I_CaL = g_CaL * m * (V - E_Ca);
    realtype I_CaT = g_CaT * x_inf(V, V_mt, s_mt) * h_T * (V - E_Ca);
    realtype I_K = g_K * n * (V - E_K);
    realtype I_SK = g_SK * s_inf(c, k_s) * (V - E_K);
    realtype I_Kir = g_Kir * x_inf(V, V_k, s_k) * (V - E_K);
    realtype I_BK = g_BK * b * (V - E_K);

    realtype v_Na_inf = x_inf(V, V_mNa, s_mNa);
    realtype v_Na_inf3 = v_Na_inf * v_Na_inf * v_Na_inf;
    realtype I_NaV = g_NaV * v_Na_inf3 * h_Na * (V - E_Na);

    realtype I_A = g_A * x_inf(V, V_a, s_a) * h * (V - E_K);
    realtype I_leak = g_leak * (V - E_leak);

    realtype I_noise = 0;//w_[0];

    realtype I = I_CaL + I_CaT + I_K + I_SK + I_Kir + I_BK + I_NaV + I_A + I_leak + I_noise;

    realtype dv = -I / C_m;

    realtype dn = (x_inf(V, V_n, s_n) - n) / tau_n;
    realtype dm = (x_inf(V, V_m, s_m) - m) / tau_m;
    realtype db = (x_inf(V, V_b, s_b) - b) / tau_BK;
    realtype dh = (x_inf(V, V_h, s_h) - h) / tau_h;
    realtype dh_T = (x_inf(V, V_ht, s_ht) - h_T) / tau_ht;
    realtype dh_Na = (x_inf(V, V_hNa, s_hNa) - h_Na) / tau_hNa;

    realtype dc = -f_c * (alpha * I_CaL + k_c * c);

    dx_[0] = dv;
    dx_[1] = dn;
    dx_[2] = dm;
    dx_[3] = db;
    dx_[4] = dh;
    dx_[5] = dh_T;
    dx_[6] = dh_Na;
    dx_[7] = dc;
    aux_[0] = 0;
}