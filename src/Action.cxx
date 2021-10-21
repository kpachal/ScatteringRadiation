#include "ScatteringRadiation/Action.h"

void Action::Build() const {

Beam * beam = new Beam();
SetUserAction(beam);

}