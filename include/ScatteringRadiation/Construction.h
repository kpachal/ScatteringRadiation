#ifndef __CONSTRUCTION__
#define __CONSTRUCTION__
#include <G4VUserDetectorConstruction.hh>
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

class G4Material;
class G4Element;

class Construction: public G4VUserDetectorConstruction
{
    public:
        Construction() {};
        virtual ~Construction() {};

        virtual G4VPhysicalVolume * Construct();
    private:

    void CreateDetector();

    G4Box * solidWorld, * solidTarget, * solidDetector;
    G4LogicalVolume * logicWorld, * logicTarget, * logicDetector;
    G4VPhysicalVolume * physWorld, * physTarget, * physDetector;   

    G4Material * vacuum, * air, * tantalum; 
};



#endif
