#include <iostream>
#include <time.h>
#include "mc.h"

using namespace std;

MonteCarlo::MonteCarlo(Parser *P)
{
	// Construction de BS
	this->mod_ = new BS(P);
	// Construction de l'option
	char *optionType;
	P->extract("option type",optionType);
	assert( (strcmp(optionType,"basket") == 0)  || (strcmp(optionType,"asian") == 0) || (strcmp(optionType,"barrier_l") == 0) ||
		(strcmp(optionType,"barrier_u") == 0) || (strcmp(optionType,"barrier") == 0) || (strcmp(optionType,"performance") == 0) );
	if ( strcmp(optionType,"basket") == 0)
		this->opt_ = new Basket(P);
	else if ( strcmp(optionType,"asian") == 0)
		this->opt_ = new Asian(P);
	else if ( strcmp(optionType,"barrier_l") == 0)
		this->opt_ = new Barrier_l(P);
	else if ( strcmp(optionType,"barrier_u") == 0)
		this->opt_ = new Barrier_u(P);
	else if ( strcmp(optionType,"barrier") == 0)
		this->opt_ = new Barrier(P);
	else
		this->opt_ = new Performance(P);
	// Initialisation du generateur a MERSENNE : type 7 page 63
	rng = pnl_rng_create(PNL_RNG_MERSENNE);
	pnl_rng_sseed(rng,time(NULL));
	// Sample number et pas de difference finie
	P->extract("sample number",this->samples_);
	this->h_ = 0.1 ;
}

MonteCarlo::~MonteCarlo()
{
#ifdef _DEBUG
	cout << "~MonteCarlo() : Ready to call pnl_rng_free ->->-> " << endl;
#endif
	pnl_rng_free(&rng);
	delete this->mod_;
	delete this->opt_;
#ifdef _DEBUG
	cout << "~MonteCarlo() : Successfull call of pnl_rng_free" << endl;
#endif
}

void MonteCarlo::price(double &prix, double &ic)
{
	PnlMat *generatedPath = pnl_mat_create(this->opt_->TimeSteps_+1,this->opt_->size_);
	double payoff = 0;
	double sumPayoff = 0;
	double sumPayoffSquare = 0;
	prix = 0;
	ic = 0;
	for(int i=0; i < this->samples_; i++)
	{
		this->mod_->asset(generatedPath, this->opt_->T_, this->opt_->TimeSteps_, this->rng);
		payoff = this->opt_->payoff(generatedPath);
		sumPayoff += payoff;
		sumPayoffSquare += payoff*payoff;
	}
	prix = exp(-this->mod_->r_*this->opt_->T_)*(sumPayoff/this->samples_);
	double x = exp(-2*this->mod_->r_*this->opt_->T_)*( (sumPayoffSquare/this->samples_) - (sumPayoff/this->samples_)*(sumPayoff/this->samples_) );
	ic = 2*1.96*x/sqrt(this->samples_);
	pnl_mat_free(&generatedPath);
}


void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic)
{
	PnlMat *generatedPath = pnl_mat_create(this->opt_->TimeSteps_+1,this->opt_->size_);
	double payoff = 0;
	double sumPayoff = 0;
	double sumPayoffSquare = 0;
	prix = 0;
	ic = 0;
	for(int i=0; i < this->samples_; i++)
	{
		this->mod_->asset(generatedPath, this->opt_->T_, this->opt_->TimeSteps_, this->rng, past);
		payoff = this->opt_->payoff(generatedPath);
		sumPayoff += payoff;
		sumPayoffSquare += payoff*payoff;
	}
	prix = exp(-this->mod_->r_*this->opt_->T_)*(sumPayoff/this->samples_);
	double x = exp(-2*this->mod_->r_*this->opt_->T_)*( (sumPayoffSquare/this->samples_) - (sumPayoff/this->samples_)*(sumPayoff/this->samples_) );
	ic = 2*1.96*x/sqrt(this->samples_);
	pnl_mat_free(&generatedPath);
}

