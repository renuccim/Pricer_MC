#include <iostream>
#include "parser.h"
#include "option.h"
#include "asian.h"
#include "basket.h"
#include "barrier_l.h"
#include "barrier_u.h"
#include "barrier.h"
#include "performance.h"
#include "pnl/pnl_matrix.h"

using namespace std;

int main(int argc, char **argv)
{
  char *infile = argv[1];
  Parser *P = new Parser(infile);
	
	char* optionType;
	P->extract("option type",optionType);
	
	Option *opt_;
	if ( strcmp(optionType,"basket") == 0)
		opt_ = new Basket(P);
	else if ( strcmp(optionType,"asian") == 0)
		opt_ = new Asian(P);
	else if ( strcmp(optionType,"barrier_l") == 0)
		opt_ = new Barrier_l(P);
	else if ( strcmp(optionType,"barrier_u") == 0)
		opt_ = new Barrier_u(P);
	else if ( strcmp(optionType,"barrier") == 0)
		opt_ = new Barrier(P);
	else
		opt_ = new Performance(P);
	opt_->print();
	
	//trajectoire test
	PnlMat *path = pnl_mat_create_from_scalar(opt_->size_,opt_->TimeSteps_+1,44);
	double payoff = opt_->payoff(path);
	
	pnl_mat_free(&path);
	delete opt_;

  exit(0);
}
