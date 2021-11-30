#include "ScatteringRadiation/Construction.h"
#include "ScatteringRadiation/SensitiveSurface.h"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

G4VPhysicalVolume * Construction::Construct()
{
	//Obtain pointer to NIST material manager
	G4NistManager* nist = G4NistManager::Instance();
	
	//Build materials for world and target
	air = nist->FindOrBuildMaterial("G4_AIR");
	vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	tantalum = nist->FindOrBuildMaterial("G4_Ta");

	//Define the world solid. The world will be 20m x 20m x20m
	G4double world_sizeX=2*CLHEP::m;
	G4double world_sizeY=2*CLHEP::m;
	G4double world_sizeZ=2*CLHEP::m;
	solidWorld = new G4Box("World",world_sizeX,world_sizeY,world_sizeZ);

	//Fill the world with air
	//G4LogicalVolume* logicWorld = 
	//	new G4LogicalVolume(solidWorld, air, "myWorld");
	// Try a vacuum world instead
	logicWorld = new G4LogicalVolume(solidWorld, vacuum, "myWorld");	

	//Create the world physical volume. The world is the only
	//physical volume with no mother volume.
	physWorld = new G4PVPlacement(0,                       //no rotation
							G4ThreeVector(),       //at (0,0,0)
							logicWorld,            //its logical volume
							"World",               //its name
							0,                     //its mother  volume
							false,                 //no boolean operation
							0,                     //copy number
							true);			       //overlaps checking                     

	//Create  the shape of a target to fire particles at
	G4double target_sizeX=2.5*CLHEP::cm;
	G4double target_sizeY=2.5*CLHEP::cm;
	G4double target_sizeZ=m_targetThickness*CLHEP::micrometer;
	solidTarget = new G4Box("Target",target_sizeX, target_sizeY, target_sizeZ);

	//Create the target logical volume by
	//assigning the material of the target to be tantalum
	logicTarget = new G4LogicalVolume(solidTarget, tantalum, "myTarget");

	//Create the target physical volume by placing it in the
	//"logicWorld" logical volume.
	physTarget = new G4PVPlacement(0,              //no rotation
							G4ThreeVector(0,0,0),  //in the center
							logicTarget,           //its logical volume
							"World",               //its name
							logicWorld,             //its mother  volume
							false,                 //no boolean operation
							0,                     //copy number
							true);			       //overlaps checking                     

    // Create the sensitive detector:
	// just panels on the outer edges of the world box
	CreateDetector();

	return physWorld;
}

void Construction::CreateDetector() {
    
	solidDetector = new G4Box("solidDetector", 2.*CLHEP::m, 2.*CLHEP::m, 0.05*CLHEP::m);

    logicDetector = new G4LogicalVolume(solidDetector, vacuum, "logicDetector");

    physDetector = new G4PVPlacement(0, G4ThreeVector(0.*CLHEP::m, 0.*CLHEP::m, 1.95*CLHEP::m), logicDetector, "physDetector", logicWorld, false, 0, true);

    //physDetector = new G4PVPlacement(0, G4ThreeVector(0.*m, 0.*m, 3.*m), logicDetector, "physDetector", logicWorld, false, 1, true);	
}

void Construction::ConstructSDandField() {

    SensitiveSurface * sens = new SensitiveSurface("SensitiveSurface");
    logicDetector->SetSensitiveDetector(sens);	
}