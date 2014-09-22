#include <iostream>
#include <time.h>
#include "bs.h"

using namespace std;

BS::BS(Parser *P)
{
	P->extract("option size",this->size_);
	P->extract("interest rate",this->r_);
	P->extract("correlation",this->rho_);
	P->extract("volatility",this->sigma_, this->size_);	
	P->extract("spot",this->spot_, this->size_);
	P->extract("trend",this->trend_, this->size_);
	// The cholesky factorization
	L = pnl_mat_create_from_scalar(this->size_,this->size_,1);
	for(int i=0; i < this->size_; i++)
	{
		for(int j=0; j < this->size_; j++)
		{
			if (i!=j)
				MLET(L,i,j) = this->rho_;
		}
	}
	pnl_mat_chol(L);
	
}

BS::~BS()
{
#ifdef _DEBUG
	cout << "~BS() : Ready to call pnl_vect_free on sigma, spot and trend ->->-> " << endl;
#endif
	pnl_vect_free(&this->sigma_);
	pnl_vect_free(&this->spot_);
	pnl_vect_free(&this->trend_);	
	pnl_mat_free(&L);
#ifdef _DEBUG
	cout << "~BS() : Successfull call of pnl_rng_free" << endl;
#endif
}

void BS::asset(PnlMat *path, double T, int N, PnlRng *rng)
{
	double step = T/N;
	double prodScal = 0;
	double sigma_d = 0;

	// The Gaussian vector
	PnlVect *G = pnl_vect_create_from_zero(this->size_);
	PnlVect *Ld = pnl_vect_create_from_zero(this->size_);	
	// First col of path contains spot
	for(int d=0; d < this->size_; d++)
		MLET(path,0,d) = GET(this->spot_,d);

	// Generation from 2 to N+1(th) column
	for(int ti=1; ti < N+1; ti++)
	{
		// Gaussian Dimension size_ generation
		pnl_vect_rng_normal(G,this->size_,rng);
		for(int d=0; d < this->size_; d++)
		{
			pnl_mat_get_row(Ld,L,d);
			prodScal = pnl_vect_scalar_prod(Ld,G);
			sigma_d = GET(this->sigma_,d);
			MLET(path,ti,d) = MGET(path,ti-1,d)*exp( (this->r_-pow(sigma_d,2)/2)*step + sigma_d*sqrt(step)*prodScal );
		}
	}
	
	// Memory free
	pnl_vect_free(&G);
	pnl_vect_free(&Ld);
}

void BS::asset(PnlMat *path, double t, int N, double T, PnlRng *rng, const PnlMat *past)
{
	int lastIndexOfPast = floor( (N/T)*t );
	double step = T/N;
	double prodScal = 0;
	double sigma_d = 0;
	// Copy of the past on the generated path
	for(int ti=0; ti < lastIndexOfPast; ti++)
	{
		for(int d=0; d < this->size_; d++)
		{
			MLET(path,ti,d) = MGET(past,ti,d);
		}
	}
	// The vector St
	PnlVect *St = pnl_vect_create_from_zero(this->size_);
	for(int d=0; d < this->size_; d++)
		pnl_vect_set(St,d,MGET(past,lastIndexOfPast,d));
	// The Gaussian vector
	PnlVect *G = pnl_vect_create_from_zero(this->size_);
	PnlVect *Ld = pnl_vect_create_from_zero(this->size_);
	// Generation from 2 to N+1(th) column
	for(int ti=lastIndexOfPast; ti < N+1; ti++)
	{
		// Gaussian Dimension size_ generation
		pnl_vect_rng_normal(G,this->size_,rng);
		for(int d=0; d < this->size_; d++)
		{
			pnl_mat_get_row(Ld,L,d);
			prodScal = pnl_vect_scalar_prod(Ld,G);
			sigma_d = GET(this->sigma_,d);
			if (ti == lastIndexOfPast)
				MLET(path,ti,d) = GET(St,d)*exp( (this->r_-pow(sigma_d,2)/2)*(ti-t) + sigma_d*sqrt(ti-t)*prodScal );
			else
				MLET(path,ti,d) = MGET(path,ti-1,d)*exp( (this->r_-pow(sigma_d,2)/2)*step + sigma_d*sqrt(step)*prodScal );
		}
	}
	// Memory free
	pnl_vect_free(&G);
	pnl_vect_free(&Ld);
	pnl_vect_free(&St);
}

void BS::asset_(PnlMat *path, double t, int N, double T, PnlRng *rng, const PnlMat *past)
{
	//pnl_mat_resize(path, N+1, size_);
	pnl_mat_set_subblock(path, past, int(N*t/T), size_);

	//pnl_mat_print(path);

	PnlVect *W = pnl_vect_create(size_);

	double Sti_1;
	double e;
	double sigma_d;
	double LdW;
	PnlVect *Ld = pnl_vect_create(size_);

	int i = int(N*t/T)+1;
	double t_  = T * (double(i)/double(N)) - t;

	pnl_vect_rng_normal(W,size_,rng);


	//cout << "before BS" << endl;
	for (int d = 0; d < size_; d++)
	{
		Sti_1 = pnl_mat_get(past, int(N*t), d);
		sigma_d = pnl_vect_get(sigma_,d);
		pnl_mat_get_row(Ld,L,d);
		LdW = pnl_vect_scalar_prod(Ld,W);
		e = exp((r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW);
		pnl_mat_set(path, i, d, Sti_1*e);
	}
	//cout << "after 1st for BS" << endl;

	t_  = T/double(N);



	for (int i = int(N*t/T)+2 ; i < N+1; i++)
	{
		// cout << ""<< endl;
		// cout << "in 2nd for BS : " << i;
		pnl_vect_rng_normal(W,size_,rng);

		for (int d = 0; d < size_; d++)
		{
			//cout << "  in 3rd for BS : " << d << endl;
			Sti_1 = pnl_mat_get(path, i-1, d);
			//cout << " Here 1 ";
			sigma_d = pnl_vect_get(sigma_,d);
			//cout << " Here 2 ";
			pnl_mat_get_row(Ld,L,d);
			//cout << " Here 3 ";
			LdW = pnl_vect_scalar_prod(Ld,W);
			//cout << " Here 4 ";
			e = exp((r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW);
			//cout << " Here 5 ";
			pnl_mat_set(path, i, d, Sti_1*e);
		}
	}

	pnl_vect_free(&Ld);
	pnl_vect_free(&W);
}

void BS::shift_asset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double T, double timestep)
{
	int pastSize = floor( (timestep/T)*t );
	for(int ti=0; ti < timestep+1; ti++)
	{
		for(int j=0; j < this->size_; j++)
		{
			if ( (ti > pastSize) && (j == d) )
				MLET(shift_path,ti,j) = (1+h)*MGET(path,ti,j);
			else
				MLET(shift_path,ti,j) = MGET(path,ti,j);
		}
	}
}

void BS::simul_market(PnlMat *path, double T, int H, PnlRng *rng)
{
	double step = T/H;
	double prodScal = 0;
	double sigma_d = 0;
	double mu_d = 0;

	// The Gaussian vector
	PnlVect *G = pnl_vect_create_from_zero(this->size_);
	PnlVect *Ld = pnl_vect_create_from_zero(this->size_);	
	// First col of path contains spot
	for(int d=0; d < this->size_; d++)
		MLET(path,0,d) = GET(this->spot_,d);

	// Generation from 2 to H+1(th) column
	for(int ti=1; ti < H+1; ti++)
	{
		// Gaussian Dimension size_ generation
		pnl_vect_rng_normal(G,this->size_,rng);
		for(int d=0; d < this->size_; d++)
		{
			pnl_mat_get_row(Ld,L,d);
			prodScal = pnl_vect_scalar_prod(Ld,G);
			sigma_d = GET(this->sigma_,d);
			mu_d = GET(this->trend_,d);
			MLET(path,ti,d) = MGET(path,ti-1,d)*exp( (mu_d-pow(sigma_d,2)/2)*step + sigma_d*sqrt(step)*prodScal );
		}
	}
	// Memory free
	pnl_vect_free(&G);
	pnl_vect_free(&Ld);
}
