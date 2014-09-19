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
		cout << " 1 - pricing à 0 " << endl;
		cout << " 2 - princing à t > 0 " << endl;
		cout << " 3 - couverture " << endl;
		cout << " 4 - exit " << endl;
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
			cout << "   2 is not implemented yet" << endl;
			cout << " ---- Pricing à t > 0 ... " << endl;
			double prix = 0;
			double ic = 0;
			//mc_->price(&prix, &ic);
			cout << " ---- prix = " << prix << endl;
			cout << " ---- ic   = " << ic << endl; 
		}else if (choix == 3){
			cout << "   3 is not implemented yet" << endl;
		}else if (choix == 4){
			break;
		}else{
			cout << "  Choix incorrect! réessayez ... " << endl;
		}
	}while(choix != 4);
	
	delete P;
	delete mc; 
  exit(0);
}
