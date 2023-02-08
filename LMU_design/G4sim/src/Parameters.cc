//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "Parameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"


//Phantom

extern const G4double worldZ1 = 100/2;                                              // half-z of world
extern const G4double worldX1 = 160/2;                                              // half-x of world
extern const G4double worldY1 = 140/2;                                              // half-y of world

extern const G4double waterZ1 = 100/2;                                              // half-z of water
extern const G4double waterX1 = 160/2;                                              // half-x of water
extern const G4double waterY1 = 140/2;                                              // half-y of water

extern const G4double phantomZ1 = 100/2;                                              // half-z of phantom
extern const G4double phantomX1 = 160/2;                                              // half-x of phantom
extern const G4double phantomY1 = 140/2;                                              // half-y of phantom
extern const G4double phantomThickness = 15/2*mm;                                     // Half wall thickness, overrides water dimensions 
                                                               
extern const G4double eWindowZ = 0.05*mm/2;                                           // half-z of entrance window
extern const G4double eWindowRadius = 7.5*mm;                                         // radius of entrance window


//Station Planes

extern const std::vector<G4double> vec1 = {0*mm, 0*mm, -worldZ1 + 47*mm};
extern const std::vector<G4double> vec2 = {0*mm, 0*mm, -worldZ1 + 48*mm};
extern const std::vector<G4double> vec3 = {0*mm, 0*mm, -worldZ1 + 49*mm};
extern const std::vector<G4double> vec4 = {0*mm, 0*mm, -worldZ1 + 50*mm};

extern const G4double theta = 90*deg;
extern const G4double cwRot = 0*deg;
extern const G4double ccRot = 90*deg;



//Fibres

extern const G4double scifiN = 33;                                                    // Number of fibres, produced 10 mm planes
extern const G4double scifiFibreRadius = 0.25*mm/2;                                   // Half-radius of fibre
extern const G4double scifiPitch = 0.300*mm;                                          // Fibre pitch, center-to-center distance
extern const G4double scifiLength = (scifiPitch*(scifiN/2)+scifiFibreRadius);         // half-length of plane
extern const G4double scifiStationDepth = scifiFibreRadius;                           // half-depth of a Station (2 Planes -> 1 Station)
extern const G4double scifiStationSide = scifiLength;                                 // Transverse face edge length of station
extern const G4double fibreLength = 10*mm;                                            // [mm] length of fibres in y direction
extern const G4double centerYPos = 0*mm;                                     // [mm]
extern const G4double centerZPos = (-50+13)*mm;                              // [mm] z-position of frame
extern const G4double sciEff = 8000;    // photons per meV
extern const G4double transEff = 0.03;  // percentage efficiency