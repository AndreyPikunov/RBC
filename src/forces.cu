#include "forces.cuh"

/******************* S T R E T C H I N G  F O R C E ***************************************/
__global__ void calculate_fl(REAL *var, size_t N, REAL *fl_x, REAL *fl_y, REAL l0, REAL kl) {

    size_t tid = get_tid();
    size_t stride = get_stride();

    REAL    l_after, l_before,
            f_after_x, f_before_x,
            f_after_y, f_before_y;


    for ( ; tid<N; tid+=stride) {

        //printf("fl : treadIdx.x = %d || blockIdx.x = %d\n", threadIdx.x, blockIdx.x);

        l_after = length(get_x(tid, var, N), get_y(tid, var, N),
                         get_x(tid+1, var, N), get_y(tid+1, var, N));
        l_before = length(get_x(tid, var, N), get_y(tid, var, N),
                          get_x(tid-1, var, N), get_y(tid-1, var, N));
        //printf("l_after = %f || l_before = %f\n",l_after , l_before);

        f_after_x = -kl*(l_after - l0)/(l_after * l0*l0) * (get_x(tid, var, N) - get_x(tid+1, var, N));
        f_before_x = -kl*(l_before - l0)/(l_before * l0*l0) * (get_x(tid, var, N) - get_x(tid-1, var, N));
        //printf("f_after_x = %f || f_before_x = %f\n",f_after_x , f_before_x);

        f_after_y = -kl*(l_after - l0)/(l_after * l0*l0) * (get_y(tid, var, N) - get_y(tid+1, var, N));
        f_before_y = -kl*(l_before - l0)/(l_before * l0*l0) * (get_y(tid, var, N) - get_y(tid-1, var, N));
        //printf("f_after_y = %f || f_before_y = %f\n",f_after_y , f_before_y);

        fl_x[tid] = f_after_x + f_before_x;
        fl_y[tid] = f_after_y + f_before_y;
        //printf("fl_x[%d] = %f || fl_y[%d] = %f\n",tid, fl_x[tid] ,tid,  fl_y[tid]);
    }
}

/******************* S W E L L I N G  F O R C E ******************************************************/
__global__ void calculate_fs(REAL *var, REAL *area, size_t N, REAL *fs_x, REAL *fs_y, REAL s0, REAL ks) {

    size_t tid = get_tid();
    size_t stride = get_stride();

    for ( ; tid<N; tid+=stride) {

        fs_x[tid] = -0.5*ks * (std::abs(*area) - s0)/(s0*s0) * (get_y(tid+1, var, N) - get_y(tid-1, var, N)) * sign(*area);
        fs_y[tid] = -0.5*ks * (std::abs(*area) - s0)/(s0*s0) * (get_x(tid-1, var, N) - get_x(tid+1, var, N)) * sign(*area);

    }
}

/******************* B E N D I N G  F O R C E ************************************/
// was copied from Hengdi code
__global__ void calculate_fb(REAL *var, size_t N, REAL *fb_x, REAL *fb_y, REAL kb) {

    size_t tid = get_tid();
    size_t stride = get_stride();

    REAL    dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4,
            kkb1, kkb2, kkb3,
            dcos1dx, dcos2dx, dcos3dx, dcos1dy, dcos2dy, dcos3dy,
            l1, l2, l1_inv, l2_inv, l3_inv, l4_inv,
            costheta1, costheta2, costheta3,
            temp;

    for ( ; tid<N; tid+=stride) {

        dx1 = get_x(tid+1, var, N) - get_x(tid, var, N);      // x1 - x
        dx2 = get_x(tid, var, N) - get_x(tid-1, var, N);      // x - x2
        dx3 = get_x(tid+2, var, N) - get_x(tid+1, var, N);    // x3 - x1
        dx4 = get_x(tid-1, var, N) - get_x(tid-2, var, N);    // x2 - x4
        dy1 = get_y(tid+1, var, N) - get_y(tid, var, N);      // y1 - y
        dy2 = get_y(tid, var, N) - get_y(tid-1, var, N);      // y - y2
        dy3 = get_y(tid+2, var, N) - get_y(tid+1, var, N);    // y3 - y1
        dy4 = get_y(tid-1, var, N) - get_y(tid-2, var, N);    // y2 - y4

        // length
        l1 = std::sqrt(dx1*dx1 + dy1*dy1);
        l2 = std::sqrt(dx2*dx2 + dy2*dy2);
        l1_inv = 1.0 / l1;
        l2_inv = 1.0 / l2;
        l3_inv = 1.0 / std::sqrt(dx3*dx3 + dy3*dy3);
        l4_inv = 1.0 / std::sqrt(dx4*dx4 + dy4*dy4);

        // trigonometrics
        costheta1 = (dx1*dx2 + dy1*dy2) * l1_inv * l2_inv;
        costheta2 = (dx3*dx1 + dy3*dy1) * l1_inv * l3_inv;
        costheta3 = (dx2*dx4 + dy2*dy4) * l2_inv * l4_inv;

        // trigonometrics derivatives
        temp = 1.0 + costheta1;
        kkb1 = kb / (temp*temp);
        temp = 1.0 + costheta2;
        kkb2 = kb / (temp*temp);
        temp = 1.0 + costheta3;
        kkb3 = kb / (temp*temp);
        REAL l13_inv = l1_inv*l1_inv*l1_inv;
        REAL l23_inv = l2_inv*l2_inv*l2_inv;
        dcos1dx =
                (dx1 - dx2) * l1_inv * l2_inv + (dx1*dx2 + dy1*dy2)*dx1 * l13_inv * l2_inv -
                (dx1*dx2 + dy1*dy2)*dx2 * l1_inv * l23_inv;
        dcos2dx = -dx3 * l3_inv * l1_inv + (dx3*dx1 + dy3*dy1)*dx1 * l3_inv * l13_inv;
        dcos3dx = dx4 * l2_inv * l4_inv - (dx2*dx4 + dy2*dy4)*dx2 * l23_inv * l4_inv;
        dcos1dy =
                (dy1 - dy2) * l1_inv * l2_inv + (dx1*dx2 + dy1*dy2)*dy1 * l13_inv * l2_inv -
                (dx1*dx2 + dy1*dy2)*dy2 * l1_inv * l23_inv;
        dcos2dy = -dy3 * l3_inv * l1_inv + (dx3*dx1 + dy3*dy1)*dy1 * l3_inv * l13_inv;
        dcos3dy = dy4 * l2_inv * l4_inv - (dx2*dx4 + dy2*dy4)*dy2 * l23_inv * l4_inv;

        // bending force
        fb_x[tid] = kkb1*dcos1dx + kkb2*dcos2dx + kkb3*dcos3dx;
        fb_y[tid] = kkb1*dcos1dy + kkb2*dcos2dy + kkb3*dcos3dy;
    }
}

//__global__ void fr(Vesicle *ves, double lr = 1.0, double kr = 1.0, Segment *segm, Vesicle *ves_out);