///*
// * Test and homeword file for the finite difference partial differential equation
// * solver classes.
// *
// */
//
//
//#include<iostream>
//#include<eigen3/Eigen/Dense>
//
///* incluces:
//#include"../linearalgebra/linearsolver.hpp"
//#include"../heatpdesolvers/heatpdesolver.hpp"
//#include"../blackscholes/gleft.hpp"
//#include"../blackscholes/gright.hpp"
//#include"../blackscholes/ftau.hpp"
//*/
//#include"../blackscholes/bspdesolver.hpp"
//#include"../blackscholes/bs.hpp"
//#include"../utils/makesolvers.hpp"
//
//using namespace std;
//using namespace Eigen;
//
//int main(){
//    cout.precision(11);
//
//    double s = 41;
//    double k = 40;
//    double vol = 0.35;
//    double t = 0.75;
//    double r = 0.04;
//    double q = 0.02;
//
//    bs closeform(s, k, vol, t, r, q);
//
//    /* European Put */
//    BsPdeSolver euro_put_fe = make_euro_put_fe(s, k, vol, t, r, q);
//    BsPdeSolver euro_put_cn_sor = make_euro_put_cn_sor(s, k, vol, t, r, q);// use default omega and tol
//    /* European Put Finite Difference Matrix and Price for alpha = 0.5, M = 4 */
//    cout << "Forward Euler European Put" << endl << euro_put_fe.solve_pde(0.5, 4) << endl
//         << "Price: " << euro_put_fe.compute_price(0.5, 4) << endl << endl;
//    cout << "Crank Nicolson Sor European Put" << endl << euro_put_cn_sor.solve_pde(0.5, 4) << endl
//         << "Price: " << euro_put_cn_sor.compute_price(0.5, 4) << endl << endl;
//
//    cout << "closed form numbers for European Put"
//         << "price: " << closeform.put_price() << endl
//         << "delta: " << closeform.put_delta() << endl
//         << "gamma: " << closeform.gamma() << endl
//         << "theta: " << closeform.put_theta() << endl << endl;
//
//
//    /* American Put */
//    BsPdeSolver amer_put_fe = make_amer_put_fe(s, k, vol, t, r, q);
//    BsPdeSolver amer_put_cn = make_amer_put_cn(s, k, vol, t, r, q);
//    /* American Put Finite Difference Matrix and Price for alpha = 0.45 and M = 4 */
//    cout << "American Put FE" << endl << amer_put_fe.solve_pde(0.45, 4) << endl
//         << "price" << endl << amer_put_fe.compute_price(0.45, 4) << endl << endl;
//    cout << "American Put CN" << endl << amer_put_cn.solve_pde(0.45, 4) << endl
//         << "price" << endl << amer_put_cn.compute_price(0.45, 4) << endl << endl;
//
//
//    /* Iterate to get the price pattern for American Put*/
//    cout << "m" << "          "
//         << "price" << "         "
//         << "delta" << "         "
//         << "gamma" << "         "
//         << "theta" << "         "
//         << endl;
//
//    cout << endl << "      American Put FE" << endl;
//    for(int m=4; m<=256; m=m*4){
//        cout << m << "  " << amer_put_fe.compute_price(0.45, m) << "  " << amer_put_fe.compute_delta(0.45, m)
//             << "  " << amer_put_fe.compute_gamma(0.45, m) << "  " << amer_put_fe.compute_theta(0.45, m) << endl;
//    }
//    cout << endl << "      American Put CN alpha 0.45" << endl;
//    for(int m=4; m<=256; m=m*4){
//        cout << m << "  " << amer_put_cn.compute_price(0.45, m) << "  " << amer_put_cn.compute_delta(0.45, m)
//             << "  " << amer_put_cn.compute_gamma(0.45, m) << "  " << amer_put_cn.compute_theta(0.45, m) << endl;
//    }
//    cout << endl << "      American Put CN alpha 5" << endl;
//    for(int m=4; m<=256; m=m*4){
//        cout << m << "  " << amer_put_cn.compute_price(5, m) << "  " << amer_put_cn.compute_delta(5, m)
//             << "  " << amer_put_cn.compute_gamma(5, m) << "  " << amer_put_cn.compute_theta(5, m) << endl;
//    }
//
//
//
//}
