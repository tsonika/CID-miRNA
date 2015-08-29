
#include "cykmatrix.h"


//The == and != operators overloaded for CykListElement keeping in view the frequent equality checks required
bool CykListElement::operator==(const CykListElement& InCle) const
{
    if (this->ruleLHS != InCle.ruleLHS) return false;
    if (this->ruleRHSpreDot != InCle.ruleRHSpreDot) return false;
    if (this->ruleRHSpostDot != InCle.ruleRHSpostDot) return false;
    if (this->score != InCle.score) return false;
    if (this->coordinates.row != InCle.coordinates.row) return false;
    if (this->coordinates.col != InCle.coordinates.col) return false;
    if (this->coordinates.type != InCle.coordinates.type) return false;
    if (this->Parents != InCle.Parents) return false;

    return true;
}

bool CykListElement::operator!=(const CykListElement& InCle) const
{
    if (*this == InCle) return false; else return true;
}

