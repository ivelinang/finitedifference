#define BOOST_TEST_MODULE mytests
#include <boost/test/included/unit_test.hpp>

#include<iostream>
#include<eigen3/Eigen/Dense>


#include "blackscholes/bspdesolver.hpp"
#include "blackscholes/bs.hpp"
#include "utils/makesolvers.hpp"

using namespace std;
using namespace Eigen;


BOOST_AUTO_TEST_CASE(myTestCase)
{
  BOOST_TEST(1 == 1);
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_CASE(pdeSolverTest)
{
	double s = 41;
	double k = 40;
	double vol = 0.35;
	double t = 0.75;
	double r = 0.04;
	double q = 0.00;

	bs closeform(s, k, vol, t, r, q);

	//European Call
	BsPdeSolver euro_call_be_lu = make_euro_call_be_lu(s, k, vol, t, r, q);

	/* European Put */
	BsPdeSolver euro_put_fe = make_euro_put_fe(s, k, vol, t, r, q);
	BsPdeSolver euro_put_cn_sor = make_euro_put_cn_sor(s, k, vol, t, r, q);// use default omega and tol

	auto Forward_Euler_European_Put = euro_put_fe.solve_pde(0.5, 4);
	auto price_fe_pde = euro_put_fe.compute_price(0.5, 4);

	auto price_be_pde = euro_call_be_lu.compute_price(0.1, 50);

	auto cn_european_put = euro_put_cn_sor.solve_pde(0.5, 4);
	auto price_cn_pde = euro_put_cn_sor.compute_price(0.5, 4);

	auto bs_price_put = closeform.put_price();
	auto bs_price_call = closeform.call_price();

		//         << "Price: " << euro_put_fe.compute_price(0.5, 4) << endl << endl;
		//    cout << "Crank Nicolson Sor European Put" << endl << euro_put_cn_sor.solve_pde(0.5, 4) << endl
		//         << "Price: " << euro_put_cn_sor.compute_price(0.5, 4) << endl << endl;
		//
		//    cout << "closed form numbers for European Put"
		//         << "price: " << closeform.put_price() << endl
		//         << "delta: " << closeform.put_delta() << endl
		//         << "gamma: " << closeform.gamma() << endl
		//         << "theta: " << closeform.put_theta() << endl << endl;

	BOOST_CHECK_CLOSE(price_be_pde, bs_price_call, 0.5);
	BOOST_CHECK_CLOSE(price_fe_pde, bs_price_put, 0.9);
	BOOST_CHECK_CLOSE(price_cn_pde, bs_price_put, 0.9);
}

