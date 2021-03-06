#include <iostream>
#include "basket.h"

using namespace std;

Basket::Basket(Parser *P)
{
	assert( P->extract("maturity", this->T_) &&
			P->extract("timestep number", this->TimeSteps_) &&
			P->extract("option size", this->size_) &&
			P->extract("option type", this->optionType_) &&
			P->extract("strike", this->K_) &&
			P->extract("payoff coefficients", this->payoffCoefficients_, this->size_) );
}

Basket::~Basket()
{
#ifdef _DEBUG
	cout << "~Basket() : Ready to call pnl_vect_free on payoff coefficients ... " << endl;
#endif
	pnl_vect_free(&this->payoffCoefficients_);
#ifdef _DEBUG
	cout << "~Basket() : Successfull call of pnl_vect_free on payoff coefficients" << endl;
#endif
}

double Basket::payoff(const PnlMat *path)
{
	PnlVect *ST = pnl_vect_create_from_zero(this->size_);
	pnl_mat_get_row(ST,path,this->TimeSteps_);
	double payoff = fmax(pnl_vect_scalar_prod(this->payoffCoefficients_,ST)-this->K_,0);
	pnl_vect_free(&ST);
	return payoff;
}

void Basket::print()
{
	std::cout << "option type " << this->optionType_ << std::endl;
	std::cout << "maturity " << this->T_ << std::endl;
	std::cout << "timestep number " << this->TimeSteps_ << std::endl;
	std::cout << "option size " << this->size_ << std::endl;
	std::cout << "strike " << this->K_ << std::endl;
	std::cout << "payoff coefficients "; pnl_vect_print_asrow(this->payoffCoefficients_);
}
