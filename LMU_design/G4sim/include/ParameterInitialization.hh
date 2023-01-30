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

#ifndef ParameterInitialization_h
#define ParameterInitialization_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

extern G4double worldZ;                                           // half-z of world
extern G4double worldX;                                           // half-x of world
extern G4double worldY;                                           // half-y of world

extern G4double waterZ;                                           // half-z of water
extern G4double waterX;                                         // half-x of water
extern G4double waterY;                                          // half-y of water

extern G4double phantomZ;                                         // half-z of phantom
extern G4double phantomX;                                        // half-x of phantom
extern G4double phantomY;                                        // half-y of phantom
extern G4double phantomThickness;                                     // Half wall thickness, overrides water dimensions 
                                                                // if togglePhantom is true
extern G4double eWindowZ;                                // half-z of entrance window
extern G4double eWindowRadius;                                 // radius of entrance window


//Station Planes
extern G4double stationN;

extern std::vector<G4double> station1Pos;                           // Station 1 depth (at station centre)
extern std::vector<G4double> station2Pos;                        // Station 2 depth (at station centre)
extern std::vector<G4double> station3Pos;                      // Station 3 depth (at station centre)
extern std::vector<G4double> station4Pos;                       // Station 4 depth (at station centre)

extern G4double theta;
extern G4double cwRot;
extern G4double ccRot;

//Fibres

extern G4double scifiN;                                                    // Number of fibres, produced 10 mm planes
extern G4double scifiFibreRadius;                                // Half-radius of fibre
extern G4double scifiPitch;                                         // Fibre pitch, center-to-center distance
extern G4double scifiLength;        // half-length of plane
extern G4double scifiStationDepth;                         // half-depth of a Station (2 Planes -> 1 Station)
extern G4double scifiStationSide;                               // Transverse face edge length of station
extern G4double fibreLength;                                            // [mm] length of fibres in y direction
extern G4double centerYPos;                                     // [mm]
extern G4double centerZPos;                             // [mm] z-position of frame

extern G4double sciEff;    // photons per meV
extern G4double transEff; // percentage efficiency
#endif
