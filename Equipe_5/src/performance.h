#ifndef _PERFORMANCE_H
#define _PERFORMANCE_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cmath>
#include "parser.h"
#include "option.h"
/// \brief Classe Option PERFORMANCE
class Performance : public Option
{
public:
  double K_; /// Strike
  PnlVect *payoffCoefficients_; /// Coefficients intervenant dans le calcul du payoff

	Performance(Parser *P);
	~Performance();

  /**
   * Calcule la valeur du payoff sur la trajectoire
   *
   * @param[in] path est une matrice de taille d x (N+1)
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * @return phi(trajectoire)
   */
  double payoff(const PnlMat *path);

	void print();
};


#endif /* _PERFORMANCE_H */
