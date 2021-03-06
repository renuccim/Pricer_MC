#include <iostream>
#include "parser.h"
#include "option.h"
#include "asian.h"
#include "basket.h"
#include "barrier_l.h"
#include "barrier_u.h"
#include "barrier.h"
#include "performance.h"
#include "bs.h"
#include "pnl/pnl_matrix.h"


using namespace std;

int main(int argc, char **argv)
{
  char *infile = argv[1];
	char *outfile = argv[2];
  Parser *P = new Parser(infile);
	
	BS *mod_ = new BS(P);
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
	// Generateur aleatoire
	PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
	pnl_rng_sseed(rng,time(NULL));


	//Simulation du passé ( on simule sur l'interval [O,T] pour pouvoir pricer à t dans [0,T] )
	PnlMat *path = pnl_mat_create_from_scalar(opt_->TimeSteps_+1,opt_->size_,0);
	mod_->asset(path,opt_->T_, opt_->TimeSteps_, rng);
	FILE *fp;
	fp = fopen(outfile, "w");
	pnl_mat_fprint(fp,path);
	fclose(fp);
	pnl_mat_free(&path);
	pnl_rng_free(&rng);
	delete opt_;
	delete mod_;
	delete P;
  exit(0);
}
