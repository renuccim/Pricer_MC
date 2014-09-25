#include <iostream>
#include "barrier.h"

using namespace std;

Barrier::Barrier(Parser *P)
{
	assert( P->extract("maturity", this->T_) &&
			P->extract("timestep number", this->TimeSteps_) &&
			P->extract("option size", this->size_) &&
			P->extract("option type", this->optionType_) &&
			P->extract("strike", this->K_) &&
			P->extract("payoff coefficients", this->payoffCoefficients_, this->size_) &&
			P->extract("lower barrier", this->lowerBarrier_, this->size_) &&
			P->extract("upper barrier", this->upperBarrier_, this->size_) );
}

Barrier::~Barrier()
{
#ifdef _DEBUG
	cout << "~Barrier() : Ready to call pnl_vect_free on payoff coefficients & lower/upper barrier  " << endl;
#endif
	pnl_vect_free(&this->lowerBarrier_);
	pnl_vect_free(&this->upperBarrier_);
	pnl_vect_free(&this->payoffCoefficients_);
#ifdef _DEBUG
	cout << "~Barrier() : Successfull call of pnl_vect_free on payoff coefficients & lower/upper barrier" << endl;
#endif
}

double Barrier::payoff(const PnlMat *path)
{
	PnlVect *ST = pnl_vect_create_from_zero(this->size_);
	bool indicatrice = true;
	for(int ti=0; ti < this->TimeSteps_+1; ti++)
	{
		pnl_mat_get_row(ST,path,ti);
		for(int d=0; d< this->size_; d++)
		{
			indicatrice = indicatrice && (GET(this->lowerBarrier_,d) <= GET(ST,d)) && (GET(this->upperBarrier_,d) >= GET(ST,d));
			if (!indicatrice)
			{
				pnl_vect_free(&ST);
				return 0;
			}
		 }
	}
	double payoff = fmax(pnl_vect_scalar_prod(this->payoffCoefficients_,ST)-this->K_ , 0);
	pnl_vect_free(&ST);
	return payoff;
}

void Barrier::print()
{
	cout << "option type" << this->optionType_ << endl;
	cout << "maturity " << this->T_ << endl;
	cout << "timestep number " << this->TimeSteps_ << endl;
	cout << "option size " << this->size_ << endl;
	cout << "K " << this->K_ << endl;
	cout << "payoff coefficients "; pnl_vect_print_asrow(this->payoffCoefficients_);
	cout << "lower barrier "; pnl_vect_print_asrow(this->lowerBarrier_);
	cout << "upper barrier "; pnl_vect_print_asrow(this->upperBarrier_);
}
