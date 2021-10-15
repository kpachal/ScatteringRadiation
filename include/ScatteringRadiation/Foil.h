#ifndef __FOIL__
#define __FOIL__
#include <G4VUserDetectorConstruction.hh>

class G4Material;
class G4Element;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

class Foil: public G4VUserDetectorConstruction
{
    public:
        Foil() {};
        virtual ~Foil() {};

        virtual G4VPhysicalVolume * Construct();
    private:
};



#endif
