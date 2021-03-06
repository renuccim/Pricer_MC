#ifndef _ASIAN_H
#define _ASIAN_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cmath>
#include "parser.h"
#include "option.h"
#include <cassert>

/// \brief Classe Option ASIAN
class Asian : public Option
{
public:
  double K_; /// Strike

	Asian(Parser *P);
	virtual ~Asian();

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


#endif /* _ASIAN_H */
