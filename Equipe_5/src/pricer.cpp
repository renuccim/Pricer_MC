#include <iostream>
#include "parser.h"
#include "mc.h"
#include "pnl/pnl_matrix.h"



using namespace std;

int main(int argc, char **argv)
{
  char *infile = argv[1];
  Parser *P = new Parser(infile);
	MonteCarlo *mc = new MonteCarlo(P);
	cout << " ------------------------------------------- " << endl;
	cout << " ------- Pricer MonteCarlo Generique ------- " << endl;
	cout << " ------------------------------------------- " << endl;
	cout << " " << endl;
	cout << " " << endl;
	int choix = 0 ;
	do{
		cout << " 1 - pricing à t = 0 " << endl;
		cout << " 2 - princing à t > 0 " << endl;
		cout << " 3 - delta à t = 0 " << endl;
		cout << " 4 - couverture " << endl;
		cout << " 5 - exit " << endl;
		cout << " Choix : " ;
		cin >> choix ;	
		if (choix == 1)
		{
			cout << " ---- Pricing à t = 0 ... " << endl;
			double prix = 0;
			double ic = 0;
			mc->price(prix, ic);
			cout << " ---- prix = " << prix << endl;
			cout << " ---- ic   = " << ic << endl; 
		}else if (choix == 2){
			cout << " ---- Pricing à t > 0 ... " << endl;
			double prix = 0;
			double ic = 0;
			cout << "      Saisir le chemin de la matrice past : exemple past/basket.past ";
			std::string chemin;
			cin >> chemin;
			cout << " " << endl;
			cout << "      Saisir t entre 0 et " << mc->opt_->T_ << " : ";
			double t = 0;
			cin >> t;
			cout << " " << endl;
			PnlMat *past = pnl_mat_create_from_file(chemin.c_str());
			mc->price(past, t, prix, ic);
			cout << " ---- prix = " << prix << endl;
			cout << " ---- ic   = " << ic << endl;
			pnl_mat_free(&past); 
		}else if (choix == 4){
			cout << " ---- Couverture ... " << endl;
			int H = 0;
			cout << "      Saisir le nombre de date de rebalancement du portefeuille de couverture H > " << mc->opt_->TimeSteps_ << " "  ;
			cin >> H;
			PnlVect *V = pnl_vect_create_from_zero(H+1);
			double PL = 0;
			mc->hedge(V,PL,H);
			pnl_vect_print_asrow(V);
			cout << " Profit & Loss : " << PL << endl;
			pnl_vect_free(&V);
		}else if (choix == 3){
			cout << "      Saisir le chemin de la matrice past : exemple past/basket.past ";
			std::string chemin;
			cin >> chemin;
			cout << " " << endl;	
			PnlMat *past = pnl_mat_create_from_file(chemin.c_str());
			PnlVect *delta = pnl_vect_create_from_zero(mc->opt_->size_);
			mc->delta(past,0,delta);
			pnl_vect_print_asrow(delta);
			pnl_vect_free(&delta);
			pnl_mat_free(&past);
		}else if (choix == 5){			
			break;
		}else{
			cout << "  Choix incorrect! réessayez ... " << endl;
		}
	}while(choix != 4);
	
	delete P;
	delete mc; 
  exit(0);
}
