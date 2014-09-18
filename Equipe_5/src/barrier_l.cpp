#include <iostream>
#include "barrier_l.h"

using namespace std;

Barrier_l::Barrier_l(Parser *P)
{
	P->extract("maturity", this->T_);
	P->extract("timestep number", this->TimeSteps_);
	P->extract("option size", this->size_);
	P->extract("option type", this->optionType_);
	P->extract("K", this->K_);
	P->extract("payoff coefficients", this->payoffCoefficients_, this->size_);
	P->extract("lower barrier", this->lowerBarrier_, this->size_);
}

Barrier_l::~Barrier_l()
{
#ifdef _DEBUG
	cout << "~Barrier_l() : Ready to call pnl_vect_free on payoff coefficients & lower barrier->->-> " << endl;
#endif
	pnl_vect_free(&this->lowerBarrier_);
	pnl_vect_free(&this->payoffCoefficients_);
#ifdef _DEBUG
	cout << "~Barrier_l() : Successfull call of pnl_vect_free on payoff coefficients & lower barrier" << endl;
#endif
}

double Barrier_l::payoff(const PnlMat *path)
{
	PnlVect *ST;
	bool indicatrice = true;
	for(int ti=1; ti <= this->TimeSteps_+1; ti++)
	{
		pnl_mat_get_col(ST,path,ti);
		for(int d=1; d<= this->size_; d++)
		{
			indicatrice = indicatrice && (GET(this->lowerBarrier_,d) <= GET(ST,d));
			if (!indicatrice)
				return 0;
		}
	}
	double payoff = fmax(pnl_vect_scalar_prod(this->payoffCoefficients_,ST)-this->K_ , 0);
	pnl_vect_free(&ST);
	return payoff;
}

void Barrier_l::print()
{
	cout << "option type" << this->optionType_ << endl;
	cout << "maturity " << this->T_ << endl;
	cout << "timestep number " << this->TimeSteps_ << endl;
	cout << "option size " << this->size_ << endl;
	cout << "K " << this->K_ << endl;
	cout << "payoff coefficients "; pnl_vect_print_asrow(this->payoffCoefficients_);
	cout << "lower barrier "; pnl_vect_print_asrow(this->lowerBarrier_);
}
