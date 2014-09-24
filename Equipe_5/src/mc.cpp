#include <iostream>
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
	cout << " ---- variance = " << x << endl;
	ic = 2*1.96*sqrt(x)/sqrt(this->samples_);
	pnl_mat_free(&generatedPath);
}


void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic)
{
	std::setprecision(9);
	PnlMat *generatedPath = pnl_mat_create(this->opt_->TimeSteps_+1,this->opt_->size_);
	double payoff = 0;
	double sumPayoff = 0;
	double sumPayoffSquare = 0;
	prix = 0;
	ic = 0;
	for(int i=0; i < this->samples_; i++)
	{
		this->mod_->asset(generatedPath, t, this->opt_->TimeSteps_, this->opt_->T_, this->rng, past);
		payoff = this->opt_->payoff(generatedPath);
		sumPayoff += payoff;
		sumPayoffSquare += payoff*payoff;
	}
	prix = exp(-this->mod_->r_*(this->opt_->T_-t))*(sumPayoff/this->samples_);
	double x = exp(-2*this->mod_->r_*(this->opt_->T_-t))*( (sumPayoffSquare/this->samples_) - (sumPayoff/this->samples_)*(sumPayoff/this->samples_) );
	cout << " ---- variance = " << x << endl;
	ic = 2*1.96*sqrt(x)/sqrt(this->samples_);
	pnl_mat_free(&generatedPath);
}

void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta)
{
	pnl_vect_set_all(delta,0);
	int pastSize = floor( (this->opt_->TimeSteps_/this->opt_->T_)*t );
	PnlMat *generatedPath = pnl_mat_create(this->opt_->TimeSteps_+1,this->opt_->size_);
	PnlMat *shiftPath_Right = pnl_mat_create(this->opt_->TimeSteps_+1,this->opt_->size_);
	PnlMat *shiftPath_Left = pnl_mat_create(this->opt_->TimeSteps_+1,this->opt_->size_);
	double tmp = 0;
	// The vector St
	PnlVect *St = pnl_vect_create_from_zero(this->opt_->size_);
	for(int d=0; d < this->opt_->size_; d++)
		pnl_vect_set(St,d,MGET(past,past->m-1,d));

	for(int j=0; j<this->samples_; j++)
	{
		//Generation d'une trajectoire
		this->mod_->asset(generatedPath, t, this->opt_->TimeSteps_, this->opt_->T_, this->rng, past);
		for(int d=0; d<this->mod_->size_; d++)
		{
			// Shift à droite
			this->mod_->shift_asset(shiftPath_Right, generatedPath, d, this->h_, t, this->opt_->T_, this->opt_->TimeSteps_);
			// Shift à droite
			this->mod_->shift_asset(shiftPath_Left, generatedPath, d, -this->h_, t, this->opt_->T_, this->opt_->TimeSteps_);
			// Delta payoff Right Left
			tmp = this->opt_->payoff(shiftPath_Right)-this->opt_->payoff(shiftPath_Left);
			pnl_vect_set(delta,d,GET(delta,d)+tmp);
		}
	}
	for(int d=0; d<this->opt_->size_; d++)
	{
		pnl_vect_set(delta,d,GET(delta,d)*exp(-this->mod_->r_*(this->opt_->T_-t))/(this->samples_*this->h_*2*GET(St,d)));
	}
	pnl_mat_free(&generatedPath);
	pnl_mat_free(&shiftPath_Right);
	pnl_mat_free(&shiftPath_Left);
	pnl_vect_free(&St);
}

void MonteCarlo::hedge(PnlVect *V, double &PL, int H, const PnlMat *marketPath)
{
	PnlVect *delta_cour = pnl_vect_create_from_zero(this->opt_->size_);
	PnlVect *delta_prec = pnl_vect_create_from_zero(this->opt_->size_);
	PnlVect *St = pnl_vect_create_from_zero(this->opt_->size_);
	double p = 0;
	double ic=0;
	double value = 0;
	
	// Affichage des valeurs du portefeuille	
	cout << " ---- \t Date n° \t PortFolio Value " << endl;
	// Initialisation du portefeuille
	this->price(p,ic);
	this->delta(marketPath,0,delta_cour);
	pnl_mat_get_row(St,marketPath,0);
	value = p-pnl_vect_scalar_prod(delta_cour,St);
	pnl_vect_set(V,0,value);
	cout << " ---- \t " << 0 << " \t " << GET(V,0) << endl; 
	// Calcul des Vi
	for(int i=1; i<H+1; i++)
	{
		pnl_vect_clone(delta_prec,delta_cour);
		this->delta(marketPath,i*this->opt_->T_/H,delta_cour);
		pnl_vect_minus_vect(delta_prec,delta_cour);
		pnl_mat_get_row(St,marketPath,i);
		value = GET(V,i-1)*exp(this->mod_->r_*this->opt_->T_/H) + pnl_vect_scalar_prod(delta_prec,St);
		pnl_vect_set(V,i,value);
		cout << " ---- \t " << i << " \t \t" << GET(V,i) << endl;
	}
	// payoff calcule le prix sur une trajectoire de taille d x (N+1)
	// => Tmp = N ; N = H ; payoff(marketPath); N = Tmp;
	int tmp = this->opt_->TimeSteps_;
	this->opt_->TimeSteps_ = H;
	double payoff = this->opt_->payoff(marketPath);
	this->opt_->TimeSteps_ = tmp;
	PL = GET(V,H) + pnl_vect_scalar_prod(delta_cour,St) - payoff;
	cout << " " << endl;
	cout << "  ---- Prix de l'option en 0 = " << p << endl;
	cout << "  ---- Erreur de couverture relative en % = " << (PL/p)*100 << endl;
	// memory free
	pnl_vect_free(&St);
	pnl_vect_free(&delta_cour);
	pnl_vect_free(&delta_prec);
}

void MonteCarlo::hedge(PnlVect *V, double &PL, int H)
{
	PnlMat *marketPath = pnl_mat_create_from_zero(H+1, this->opt_->size_);
	this->mod_->asset(marketPath, this->opt_->T_,H , this->rng);
	this->hedge(V,PL,H,marketPath);
	pnl_mat_free(&marketPath);
}
