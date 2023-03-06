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

#ifndef Parameters_h
#define Parameters_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

extern const G4double worldZ1;                                           // half-z of world
extern const G4double worldX1;                                           // half-x of world
extern const G4double worldY1;                                           // half-y of world

extern const G4double waterZ1;                                           // half-z of water
extern const G4double waterX1;                                         // half-x of water
extern const G4double waterY1;                                          // half-y of water

extern const G4double phantomZ1;                                         // half-z of phantom
extern const G4double phantomX1;                                        // half-x of phantom
extern const G4double phantomY1;                                        // half-y of phantom
extern const G4double phantomThickness1;                                     // Half wall thickness, overrides water dimensions 
                                                                // if togglePhantom is true
extern const G4double eWindowZ1;                                // half-z of entrance window
extern const G4double eWindowRadius1;                                 // radius of entrance window


//Station Planes
extern const G4double theta1;
extern const G4double cwRot1;
extern const G4double ccRot1;

extern const std::vector<double> vec1;
extern const std::vector<double> vec2;
extern const std::vector<double> vec3;
extern const std::vector<double> vec4;

//Fibres

extern const G4double scifiN;                                                    // Number of fibres, produced 10 mm planes
extern const G4double scifiFibreRadius;                                // Half-radius of fibre
extern const G4double scifiPitch;                                         // Fibre pitch, center-to-center distance
extern const G4double scifiLength;        // half-length of plane
extern const G4double scifiStationDepth;                         // half-depth of a Station (2 Planes -> 1 Station)
extern const G4double scifiStationSide;                               // Transverse face edge length of station
extern const G4double fibreLength;                                            // [mm] length of fibres in y direction
extern const G4double centerYPos;                                     // [mm]
extern const G4double centerZPos;                             // [mm] z-position of frame

extern const G4double sciEff;    // photons per meV
extern const G4double transEff; // percentage efficiency
#endif


    
