#include <iostream>
#include "performance.h"

using namespace std;

Performance::Performance(Parser *P)
{
	assert( P->extract("maturity", this->T_) &&
			P->extract("timestep number", this->TimeSteps_) &&
			P->extract("option size", this->size_) &&
			P->extract("option type", this->optionType_) &&
			P->extract("payoff coefficients", this->payoffCoefficients_, this->size_) );
}

Performance::~Performance()
{
#ifdef _DEBUG
	cout << "~Performance() : Ready to call pnl_vect_free on payoff coefficients " << endl;
#endif
	pnl_vect_free(&this->payoffCoefficients_);
#ifdef _DEBUG
	cout << "~Performance() : Successfull call of pnl_vect_free on payoff coefficients " << endl;
#endif
}

double Performance::payoff(const PnlMat *path)
{
	PnlVect *STi = pnl_vect_create_from_zero(this->size_);
	PnlVect *STi_1 = pnl_vect_create_from_zero(this->size_);
	double payoff = 0;
	for(int i=1; i <= this->TimeSteps_; i++)
	{
		pnl_mat_get_row(STi,path,i);
		pnl_mat_get_row(STi_1,path,i-1);
		payoff += pnl_vect_scalar_prod(this->payoffCoefficients_,STi)/pnl_vect_scalar_prod(this->payoffCoefficients_,STi_1);
	}
	payoff = payoff/this->TimeSteps_ - 1;
	payoff = 1 + fmin(payoff,0.1);
	pnl_vect_free(&STi);
	pnl_vect_free(&STi_1);
	return payoff;
}

void Performance::print()
{
	cout << "option type" << this->optionType_ << endl;
	cout << "maturity " << this->T_ << endl;
	cout << "timestep number " << this->TimeSteps_ << endl;
	cout << "option size " << this->size_ << endl;
	cout << "payoff coefficients "; pnl_vect_print_asrow(this->payoffCoefficients_);
}
