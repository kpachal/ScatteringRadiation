#ifndef __CONSTRUCTION__
#define __CONSTRUCTION__
#include <G4VUserDetectorConstruction.hh>
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

enum target_options {foil, wire, gas};

class G4Material;
class G4Element;

class Construction: public G4VUserDetectorConstruction
{
public:
        Construction(double target_thickness, target_options target_type) { m_targetThickness = target_thickness; m_target_type = target_type; };
        virtual ~Construction() {};

        virtual G4VPhysicalVolume * Construct();
private:

    void CreateDetector();

    //G4Box * solidWorld, * solidTarget, * solidDetector;
    G4VSolid * solidWorld, * solidTarget, * solidDetector;
    G4LogicalVolume * logicWorld, * logicTarget, * logicDetector;
    G4VPhysicalVolume * physWorld, * physTarget, * physDetector;   

    G4Material * vacuum, * air, * tantalum; 
    virtual void ConstructSDandField();

    double m_targetThickness;
    target_options m_target_type;

};



#endif
