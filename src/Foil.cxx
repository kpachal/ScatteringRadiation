#include "ScatteringRadiation/Foil.h"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"

G4VPhysicalVolume * Foil::Construct()
{
	//Obtain pointer to NIST material manager
	G4NistManager* nist = G4NistManager::Instance();
	
	//Build materials for world and target
	G4Material* Pb =  nist->FindOrBuildMaterial("G4_Pb");
	G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* Ta = nist->FindOrBuildMaterial("G4_Ta");

	//Define the world solid. The world will be 20m x 20m x20m
	G4double world_sizeX=2*CLHEP::m;
	G4double world_sizeY=2*CLHEP::m;
	G4double world_sizeZ=2*CLHEP::m;
	G4Box* solidWorld = 
		new G4Box("World",world_sizeX,world_sizeY,world_sizeZ);

	//Fill the world with air
	G4LogicalVolume* logicWorld = 
		new G4LogicalVolume(solidWorld, air, "myWorld");

	//Create the world physical volume. The world is the only
	//physical volume with no mother volume.
	G4VPhysicalVolume* physWorld = 
		new G4PVPlacement(0,                       //no rotation
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
	G4double target_sizeZ=1*CLHEP::micrometer;
	G4Box* solidTarget = 
		new G4Box("Target",target_sizeX, target_sizeY, target_sizeZ);

	//Create the target logical volume by
	//assigning the material of the target to be Pb
//	G4LogicalVolume* logicTarget = 
//		new G4LogicalVolume(solidTarget, Pb, "myTarget");
	G4LogicalVolume* logicTarget = 
		new G4LogicalVolume(solidTarget, Ta, "myTarget");

	//Create the target physical volume by placing it in the
	//"logicWorld" logical volume.
	G4VPhysicalVolume* physTarget = 
		new G4PVPlacement(0,                       //no rotation
							G4ThreeVector(),       //at (0,0,0)
							logicTarget,           //its logical volume
							"World",               //its name
							logicWorld,             //its mother  volume
							false,                 //no boolean operation
							0,                     //copy number
							true);			       //overlaps checking                     

	return physWorld;
}
