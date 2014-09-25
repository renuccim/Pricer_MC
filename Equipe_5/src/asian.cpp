#include <iostream>
#include "asian.h"

using namespace std;

Asian::Asian(Parser *P)
{
	assert( P->extract("maturity", this->T_) &&
			P->extract("timestep number", this->TimeSteps_) &&
			P->extract("option size", this->size_) &&
			P->extract("option type", this->optionType_) &&
			P->extract("strike", this->K_) );
}

Asian::~Asian()
{
#ifdef _DEBUG
	cout << "~Asian()" << endl;
#endif
}

double Asian::payoff(const PnlMat *path)
{
	// Dimension D = 1
	double payoff = fmax( (pnl_mat_sum(path)/(this->TimeSteps_+1)) - this->K_ , 0);
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
