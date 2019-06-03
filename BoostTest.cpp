#define BOOST_TEST_MODULE mytests
#include <boost/test/included/unit_test.hpp>

#include<iostream>
#include<eigen3/Eigen/Dense>


#include "blackscholes/bspdesolver.hpp"
#include "blackscholes/bs.hpp"
#include "utils/makesolvers.hpp"

#include "payoff.h"
#include "vanillaoption.h"
#include "pde.h"
#include "fdm.h"

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


BOOST_AUTO_TEST_CASE(pdeSolverBetterTest)
{

	double K = 0.5;  // Strike price
	double r = 0.05;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 1.00;    // One year until expiry

						// FDM discretisation parameters
	double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
	unsigned long J = 20;
	double t_dom = T;         // Time period as for the option
	unsigned long N = 20;

	// Create the PayOff and Option objects
	PayOffCall pay_off_call = PayOffCall(K);
	VanillaOption call_option = VanillaOption(K, r, T, v, &pay_off_call);

	// Create the PDE and FDM objects
	BlackScholesPDE bs_pde = BlackScholesPDE(&call_option);
	FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, &bs_pde);

	// Run the FDM solver
	fdm_euler.step_march();

	// Delete the PDE, PayOff and Option objects
	//delete bs_pde;
	//delete call_option;
	//delete pay_off_call;

	BOOST_CHECK(TRUE);

}

BOOST_AUTO_TEST_CASE(pdeSolverBetterTest_B) {
	double S = 0.40;
	double K = 0.40;
	double v = 0.35;
	double T = 1.0;
	double r = 0.04;
	double q = 0.00;

	bs closeform(S, K, v, T, r, q);


	// FDM discretisation parameters
	double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
	unsigned long J = 20;
	double t_dom = T;         // Time period as for the option
	unsigned long N = 50;

	// Create the PayOff and Option objects
	PayOffCall pay_off_call = PayOffCall(K);
	VanillaOption call_option = VanillaOption(K, r, T, v, &pay_off_call);

	// Create the PDE and FDM objects
	BlackScholesPDE bs_pde = BlackScholesPDE(&call_option);
	FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, &bs_pde);

	//fdm_euler.step_march();
	double price_pde = fdm_euler.calcPrice(S);
	double price_bs = closeform.call_price();

	// Delete the PDE, PayOff and Option objects
	//delete bs_pde;
	//delete call_option;
	//delete pay_off_call;

	BOOST_CHECK_CLOSE(price_pde, price_bs, 0.5);
}

