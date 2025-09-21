// MolecularSpectroscopy.cc
// Geant4 simulation for molecular light spectroscopy and scattering
// Compile with: g++ -o spectroscopy MolecularSpectroscopy.cc `geant4-config --cflags --libs`

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4OpticalPhysics.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4ParticleGun.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include <vector>
#include <fstream>
#include <map>
#include <cmath>

// ============= Spectral Data Storage =============
class SpectralData {
public:
    std::map<G4double, G4int> absorptionSpectrum;
    std::map<G4double, G4int> emissionSpectrum;
    std::map<G4double, G4int> scatteringSpectrum;
    std::vector<G4double> rayleighEvents;
    
    void RecordAbsorption(G4double wavelength) {
        absorptionSpectrum[wavelength]++;
    }
    
    void RecordEmission(G4double wavelength) {
        emissionSpectrum[wavelength]++;
    }
    
    void RecordScattering(G4double wavelength) {
        scatteringSpectrum[wavelength]++;
    }
    
    void SaveToFile(const G4String& filename) {
        std::ofstream file(filename);
        
        file << "# Molecular Spectroscopy Simulation Results\n";
        file << "# Wavelength(nm), Absorption, Emission, Scattering\n";
        
        std::set<G4double> allWavelengths;
        for(auto& p : absorptionSpectrum) allWavelengths.insert(p.first);
        for(auto& p : emissionSpectrum) allWavelengths.insert(p.first);
        for(auto& p : scatteringSpectrum) allWavelengths.insert(p.first);
        
        for(auto wl : allWavelengths) {
            file << wl/nm << " " 
                 << absorptionSpectrum[wl] << " "
                 << emissionSpectrum[wl] << " "
                 << scatteringSpectrum[wl] << "\n";
        }
        
        file.close();
        G4cout << "Spectral data saved to " << filename << G4endl;
    }
};

static SpectralData* spectralData = new SpectralData();

// ============= Molecular Material Definition =============
class MolecularMaterial {
public:
    static G4Material* CreateMolecularSample(const G4String& name, 
                                            G4double density,
                                            G4double absorptionLength = 1.0*cm,
                                            G4double rayleighLength = 10.0*cm) {
        
        G4NistManager* nist = G4NistManager::Instance();
        
        // Base material (e.g., organic molecule in solution)
        G4Material* baseMat = new G4Material(name, density, 3);
        baseMat->AddElement(nist->FindOrBuildElement("C"), 6);
        baseMat->AddElement(nist->FindOrBuildElement("H"), 12);
        baseMat->AddElement(nist->FindOrBuildElement("O"), 6);
        
        // Define optical properties
        const G4int nEntries = 100;
        G4double photonEnergy[nEntries];
        G4double refractiveIndex[nEntries];
        G4double absorption[nEntries];
        G4double rayleigh[nEntries];
        G4double wlsAbsorption[nEntries];
        G4double wlsEmission[nEntries];
        
        // Generate wavelength-dependent properties (350-750 nm range)
        for(G4int i = 0; i < nEntries; i++) {
            G4double wavelength = (350.0 + i * 4.0) * nm;
            photonEnergy[i] = h_Planck * c_light / wavelength;
            
            // Refractive index (Cauchy equation approximation)
            refractiveIndex[i] = 1.5 + 0.01/(wavelength/nm/500.0) + 
                                 0.001/pow(wavelength/nm/500.0, 2);
            
            // Absorption spectrum (Gaussian peak around 450nm)
            G4double absCenter = 450.0 * nm;
            G4double absWidth = 30.0 * nm;
            absorption[i] = absorptionLength * 
                           exp(-pow((wavelength - absCenter)/absWidth, 2));
            
            // Rayleigh scattering (λ^-4 dependence)
            rayleigh[i] = rayleighLength * pow(wavelength/400.0/nm, 4);
            
            // Fluorescence absorption (peak at 450nm)
            wlsAbsorption[i] = 100.0*m * 
                              exp(-pow((wavelength - 450.0*nm)/(30.0*nm), 2));
            
            // Fluorescence emission (Stokes shifted to 520nm)
            wlsEmission[i] = exp(-pow((wavelength - 520.0*nm)/(40.0*nm), 2));
        }
        
        G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);
        mpt->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries);
        mpt->AddProperty("RAYLEIGH", photonEnergy, rayleigh, nEntries);
        mpt->AddProperty("WLSABSLENGTH", photonEnergy, wlsAbsorption, nEntries);
        mpt->AddProperty("WLSCOMPONENT", photonEnergy, wlsEmission, nEntries);
        
        // Fluorescence quantum yield and time constant
        mpt->AddConstProperty("WLSTIMECONSTANT", 8.0*ns);
        mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.8); // Quantum yield
        
        baseMat->SetMaterialPropertiesTable(mpt);
        
        return baseMat;
    }
    
    static G4Material* CreateSolvent() {
        G4NistManager* nist = G4NistManager::Instance();
        G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
        
        const G4int nEntries = 50;
        G4double photonEnergy[nEntries];
        G4double refractiveIndex[nEntries];
        G4double absorption[nEntries];
        
        for(G4int i = 0; i < nEntries; i++) {
            G4double wavelength = (300.0 + i * 10.0) * nm;
            photonEnergy[i] = h_Planck * c_light / wavelength;
            refractiveIndex[i] = 1.333;
            absorption[i] = 100.0*m; // Very low absorption
        }
        
        G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);
        mpt->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries);
        
        water->SetMaterialPropertiesTable(mpt);
        return water;
    }
};

// ============= Detector Construction =============
class SpectroscopyDetectorConstruction : public G4VUserDetectorConstruction {
public:
    virtual G4VPhysicalVolume* Construct() {
        // World volume
        G4NistManager* nist = G4NistManager::Instance();
        G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
        
        G4Box* worldBox = new G4Box("World", 50*cm, 50*cm, 50*cm);
        G4LogicalVolume* worldLV = new G4LogicalVolume(worldBox, air, "World");
        G4VPhysicalVolume* worldPV = new G4PVPlacement(0, G4ThreeVector(), 
                                                       worldLV, "World", 
                                                       0, false, 0);
        
        // Sample cuvette (quartz)
        G4Material* quartz = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
        AddQuartzProperties(quartz);
        
        G4Box* cuvetteOuter = new G4Box("CuvetteOuter", 1.25*cm, 2.5*cm, 2.5*cm);
        G4Box* cuvetteInner = new G4Box("CuvetteInner", 1.0*cm, 2.3*cm, 2.3*cm);
        
        G4SubtractionSolid* cuvetteSolid = 
            new G4SubtractionSolid("Cuvette", cuvetteOuter, cuvetteInner);
        
        G4LogicalVolume* cuvetteLV = 
            new G4LogicalVolume(cuvetteSolid, quartz, "Cuvette");
        
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), 
                         cuvetteLV, "Cuvette", worldLV, false, 0);
        
        // Sample solution
        G4Material* sample = MolecularMaterial::CreateMolecularSample(
            "OrganicMolecule", 1.2*g/cm3, 0.5*cm, 5.0*cm);
        
        G4LogicalVolume* sampleLV = 
            new G4LogicalVolume(cuvetteInner, sample, "Sample");
        
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                         sampleLV, "Sample", worldLV, false, 0);
        
        // Detector array (photodiode simulation)
        G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
        
        // Main detector (transmitted light)
        G4Box* detectorBox = new G4Box("Detector", 2*cm, 2*cm, 0.5*mm);
        G4LogicalVolume* detectorLV = 
            new G4LogicalVolume(detectorBox, silicon, "Detector");
        
        new G4PVPlacement(0, G4ThreeVector(5*cm, 0, 0),
                         detectorLV, "TransmissionDetector", worldLV, false, 0);
        
        // 90-degree detector (scattering)
        new G4PVPlacement(0, G4ThreeVector(0, 5*cm, 0),
                         detectorLV, "ScatteringDetector", worldLV, false, 1);
        
        // Fluorescence detector (opposite side)
        new G4PVPlacement(0, G4ThreeVector(-5*cm, 0, 0),
                         detectorLV, "FluorescenceDetector", worldLV, false, 2);
        
        return worldPV;
    }
    
private:
    void AddQuartzProperties(G4Material* quartz) {
        const G4int nEntries = 50;
        G4double photonEnergy[nEntries];
        G4double refractiveIndex[nEntries];
        G4double absorption[nEntries];
        
        for(G4int i = 0; i < nEntries; i++) {
            G4double wavelength = (200.0 + i * 12.0) * nm;
            photonEnergy[i] = h_Planck * c_light / wavelength;
            refractiveIndex[i] = 1.46;
            absorption[i] = 50.0*m; // Very transparent
        }
        
        G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);
        mpt->AddProperty("ABSLENGTH", photonEnergy, absorption, nEntries);
        
        quartz->SetMaterialPropertiesTable(mpt);
    }
};

// ============= Primary Generator (Light Source) =============
class SpectroscopyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
    G4ParticleGun* fParticleGun;
    G4bool fBroadband;
    
public:
    SpectroscopyPrimaryGeneratorAction(G4bool broadband = false) 
        : fBroadband(broadband) {
        fParticleGun = new G4ParticleGun(1);
        fParticleGun->SetParticleDefinition(G4OpticalPhoton::Definition());
        fParticleGun->SetParticlePosition(G4ThreeVector(-10*cm, 0, 0));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1, 0, 0));
        fParticleGun->SetParticlePolarization(G4ThreeVector(0, 1, 0));
    }
    
    virtual ~SpectroscopyPrimaryGeneratorAction() {
        delete fParticleGun;
    }
    
    virtual void GeneratePrimaries(G4Event* event) {
        G4double wavelength;
        
        if(fBroadband) {
            // Generate broadband spectrum (350-750 nm)
            wavelength = (350.0 + G4UniformRand() * 400.0) * nm;
        } else {
            // Monochromatic source (e.g., laser at 488 nm)
            wavelength = 488.0 * nm;
        }
        
        G4double energy = h_Planck * c_light / wavelength;
        fParticleGun->SetParticleEnergy(energy);
        
        // Add slight beam divergence
        G4double theta = G4RandGauss::shoot(0, 0.01);
        G4double phi = G4UniformRand() * 2.0 * pi;
        G4ThreeVector direction(cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi));
        fParticleGun->SetParticleMomentumDirection(direction);
        
        fParticleGun->GeneratePrimaryVertex(event);
    }
};

// ============= Stepping Action (Process Tracking) =============
class SpectroscopySteppingAction : public G4UserSteppingAction {
public:
    virtual void UserSteppingAction(const G4Step* step) {
        G4Track* track = step->GetTrack();
        
        if(track->GetDefinition() != G4OpticalPhoton::Definition()) return;
        
        const G4VProcess* process = 
            step->GetPostStepPoint()->GetProcessDefinedStep();
        
        if(!process) return;
        
        G4String processName = process->GetProcessName();
        G4double wavelength = h_Planck * c_light / track->GetKineticEnergy();
        
        if(processName == "OpAbsorption") {
            spectralData->RecordAbsorption(wavelength);
        } else if(processName == "OpRayleigh") {
            spectralData->RecordScattering(wavelength);
            spectralData->rayleighEvents.push_back(wavelength);
        } else if(processName == "OpWLS") {
            // Record both absorption and emission
            spectralData->RecordAbsorption(wavelength);
            // Get emitted wavelength from next step
            if(track->GetTrackStatus() != fStopAndKill) {
                G4double emittedWL = h_Planck * c_light / 
                                    track->GetKineticEnergy();
                spectralData->RecordEmission(emittedWL);
            }
        }
        
        // Check if photon reached detector
        G4String volumeName = 
            step->GetPostStepPoint()->GetTouchableHandle()
            ->GetVolume()->GetName();
        
        if(volumeName.contains("Detector")) {
            // Record detected photon
            if(volumeName == "TransmissionDetector") {
                // Direct transmission
            } else if(volumeName == "ScatteringDetector") {
                spectralData->RecordScattering(wavelength);
            } else if(volumeName == "FluorescenceDetector") {
                spectralData->RecordEmission(wavelength);
            }
        }
    }
};

// ============= Action Initialization =============
class SpectroscopyActionInitialization : public G4VUserActionInitialization {
private:
    G4bool fBroadband;
    
public:
    SpectroscopyActionInitialization(G4bool broadband = false) 
        : fBroadband(broadband) {}
    
    virtual void BuildForMaster() const {}
    
    virtual void Build() const {
        SetUserAction(new SpectroscopyPrimaryGeneratorAction(fBroadband));
        SetUserAction(new SpectroscopySteppingAction());
    }
};

// ============= Physics List =============
class SpectroscopyPhysicsList : public FTFP_BERT {
public:
    SpectroscopyPhysicsList() : FTFP_BERT() {
        // Add optical physics
        G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
        opticalPhysics->SetWLSTimeProfile("exponential");
        opticalPhysics->SetScintillationYieldFactor(1.0);
        opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
        opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);
        
        // Configure processes
        opticalPhysics->Configure(kWLS, true);
        opticalPhysics->Configure(kCerenkov, false); // No Cerenkov for this app
        opticalPhysics->Configure(kScintillation, true);
        opticalPhysics->Configure(kAbsorption, true);
        opticalPhysics->Configure(kRayleigh, true);
        opticalPhysics->Configure(kMieHG, false);
        opticalPhysics->Configure(kBoundary, true);
        
        RegisterPhysics(opticalPhysics);
        
        // High precision EM physics
        ReplacePhysics(new G4EmStandardPhysics_option4());
    }
};

// ============= Main Function =============
int main(int argc, char** argv) {
    // Parse command line arguments
    G4bool broadband = false;
    G4int nEvents = 10000;
    G4String outputFile = "spectroscopy_results.txt";
    
    for(int i = 1; i < argc; i++) {
        G4String arg = argv[i];
        if(arg == "--broadband") {
            broadband = true;
        } else if(arg == "--events" && i+1 < argc) {
            nEvents = std::atoi(argv[++i]);
        } else if(arg == "--output" && i+1 < argc) {
            outputFile = argv[++i];
        }
    }
    
    G4cout << "==================================" << G4endl;
    G4cout << "Molecular Spectroscopy Simulation" << G4endl;
    G4cout << "==================================" << G4endl;
    G4cout << "Light source: " << (broadband ? "Broadband" : "Monochromatic") << G4endl;
    G4cout << "Number of events: " << nEvents << G4endl;
    G4cout << "Output file: " << outputFile << G4endl;
    
    // Construct the run manager
    G4RunManager* runManager = new G4RunManager;
    
    // Set mandatory initialization classes
    runManager->SetUserInitialization(new SpectroscopyDetectorConstruction());
    runManager->SetUserInitialization(new SpectroscopyPhysicsList());
    runManager->SetUserInitialization(new SpectroscopyActionInitialization(broadband));
    
    // Initialize G4 kernel
    runManager->Initialize();
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if(argc == 1) {
        // Interactive mode with visualization
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
        
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        UImanager->ApplyCommand("/vis/open OGL 800x600-0+0");
        UImanager->ApplyCommand("/vis/drawVolume");
        UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 70 20");
        UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
        
        ui->SessionStart();
        
        delete ui;
        delete visManager;
    } else {
        // Batch mode
        G4String command = "/run/beamOn " + std::to_string(nEvents);
        UImanager->ApplyCommand(command);
    }
    
    // Save results
    spectralData->SaveToFile(outputFile);
    
    // Cleanup
    delete spectralData;
    delete runManager;
    
    return 0;
}

/* 
Compilation Instructions:
-------------------------
1. Ensure Geant4 is installed and configured
2. Create a CMakeLists.txt file:

cmake_minimum_required(VERSION 3.16)
project(MolecularSpectroscopy)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

add_executable(spectroscopy MolecularSpectroscopy.cc)
target_link_libraries(spectroscopy ${Geant4_LIBRARIES})

3. Build:
   mkdir build && cd build
   cmake ..
   make

4. Run:
   ./spectroscopy --broadband --events 100000 --output spectrum.txt

Features:
---------
- Molecular absorption and emission spectra
- Rayleigh scattering simulation
- Fluorescence (wavelength shifting)
- Multiple detector configuration
- Broadband and monochromatic light sources
- Wavelength-dependent material properties
- Spectral data recording and analysis

Output:
-------
The simulation generates a text file with wavelength-resolved:
- Absorption spectrum
- Emission spectrum (fluorescence)
- Scattering spectrum
- Individual photon tracking data

This data can be imported into analysis software for:
- Spectral peak identification
- Quantum yield calculation
- Scattering cross-section determination
- Molecular fingerprinting
*/