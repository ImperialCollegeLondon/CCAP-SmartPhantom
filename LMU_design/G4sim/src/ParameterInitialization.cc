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

#include "ParameterInitialization.hh"
#include "G4SystemOfUnits.hh"

//Phantom

G4double worldZ = 100*mm/2;                                              // half-z of world
G4double worldX = 160*mm/2;                                              // half-x of world
G4double worldY = 140*mm/2;                                              // half-y of world

G4double waterZ = 100*mm/2;                                              // half-z of water
G4double waterX = 160*mm/2;                                              // half-x of water
G4double waterY = 140*mm/2;                                              // half-y of water

G4double phantomZ = waterZ;                                              // half-z of phantom
G4double phantomX = waterX;                                              // half-x of phantom
G4double phantomY = waterY;                                              // half-y of phantom
G4double phantomThickness = 15/2*mm;                                     // Half wall thickness, overrides water dimensions 
                                                                // if togglePhantom is true
G4double eWindowZ = 0.05*mm/2;                                           // half-z of entrance window
G4double eWindowRadius = 7.5*mm;                                         // radius of entrance window


//Station Planes

std::vector<G4double> vec1 = {0*mm, 0*mm, -waterZ + 20*mm};
std::vector<G4double> vec2 = {0*mm, 0*mm, -waterZ + 22*mm};
std::vector<G4double> vec3 = {0*mm, 0*mm, -waterZ + 24*mm};
std::vector<G4double> vec4 = {0*mm, 0*mm, -waterZ + 26*mm};

std::vector<G4double> station1Pos = vec1;                                             // Station 1 depth (at station centre)
std::vector<G4double> station2Pos = vec2;                                             // Station 2 depth (at station centre)
std::vector<G4double> station3Pos = vec3;                                             // Station 3 depth (at station centre)
std::vector<G4double> station4Pos = vec4;                                             // Station 4 depth (at station centre)

G4double theta = 90*deg;
G4double cwRot = 0*deg;
G4double ccRot = 90*deg;

//Fibres

G4double scifiN = 33;                                                    // Number of fibres, produced 10 mm planes
G4double scifiFibreRadius = 0.25*mm/2;                                   // Half-radius of fibre
G4double scifiPitch = 0.300*mm;                                          // Fibre pitch, center-to-center distance
G4double scifiLength = (scifiPitch*(scifiN/2)+scifiFibreRadius);         // half-length of plane
G4double scifiStationDepth = scifiFibreRadius;                           // half-depth of a Station (2 Planes -> 1 Station)
G4double scifiStationSide = scifiLength;                                 // Transverse face edge length of station
G4double fibreLength = 10*mm;                                            // [mm] length of fibres in y direction
G4double centerYPos = 0*mm;                                     // [mm]
G4double centerZPos = (-50+13)*mm;                              // [mm] z-position of frame
G4double sciEff = 8000;    // photons per meV
G4double transEff = 0.03;  // percentage efficiency