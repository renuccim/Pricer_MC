#include <iostream>
#include "asian.h"

using namespace std;

Asian::Asian(Parser *P)
{
	P->extract("maturity", this->T_);
	P->extract("timestep number", this->TimeSteps_);
	P->extract("option size", this->size_);
	P->extract("option type", this->optionType_);
	P->extract("K", this->K_);
}

Asian::~Asian()
{
#ifdef _DEBUG
	cout << "~Asian()" << endl;
#endif
}

double Asian::payoff(const PnlMat *path)
{
	PnlVect *ST;
	// Dimension D = 1
	pnl_mat_get_row(ST,path,1);
	PnlVect *Id = pnl_vect_create_from_scalar(this->TimeSteps_+1,1);
	double payoff = fmax( (pnl_vect_scalar_prod(Id,ST)/this->TimeSteps_) - this->K_ , 0);
	pnl_vect_free(&ST);
	pnl_vect_free(&Id);
	return payoff;
}

void Asian::print()
{
	cout << "option type" << this->optionType_ << endl;
	cout << "maturity " << this->T_ << endl;
	cout << "timestep number " << this->TimeSteps_ << endl;
	cout << "option size " << this->size_ << endl;
	cout << "K " << this->K_ << endl;
}