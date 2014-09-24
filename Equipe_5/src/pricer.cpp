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
		cout << " 2 - pricing à t > 0 " << endl;
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
			cout << " ---- ic   = " << ic << endl; 
			cout << " ---- prix = " << prix << endl;
		}else if (choix == 2){
			cout << " ---- Pricing à t > 0 ... " << endl;
			double prix = 0;
			double ic = 0;
			cout << "      Saisir le chemin de la matrice past : exemple past/basket.past ";
			std::string chemin;
			cin >> chemin;
			cout << " " << endl;
			cout << "      Saisir t dans [0, " << mc->opt_->T_ << "] : ";
			double t = 0;
			cin >> t;
			cout << " " << endl;
			if (t > mc->opt_->T_)
			{
				cout << "      Erreur : t = " << t << " > Maturité de l'option !" << endl;
				continue;               
			}
			PnlMat *past = pnl_mat_create_from_file(chemin.c_str());
			mc->price(past, t, prix, ic);
			cout << " ---- ic   = " << ic << endl;
			cout << " ---- prix = " << prix << endl;
			pnl_mat_free(&past); 
		}else if (choix == 4){
			cout << " ---- Couverture ... " << endl;
			int H = 0;
			cout << "      Saisir le nombre de date de rebalancement du portefeuille de couverture H qui soit multiple de " << mc->opt_->TimeSteps_ << endl  ;
			cout << "      Typiquement vous pouvez choisir une date par jour ou par semaine . H : " ;
			cin >> H;
			if (H < mc->opt_->TimeSteps_)
			{
				cout << "      Warning : H = " << H << " < TimeSteps " << endl;
				continue;
			}
			PnlVect *V = pnl_vect_create_from_zero(H+1);
			double PL = 0;
			mc->hedge(V,PL,H);
			cout << " Profit & Loss : " << PL << endl;
			pnl_vect_free(&V);
		}else if (choix == 3){
			PnlMat *past = pnl_mat_create(1,mc->opt_->size_);
			pnl_mat_set_row(past,mc->mod_->spot_,0);
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
